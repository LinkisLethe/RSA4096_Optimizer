/*****************************************************************************
Filename    : rsa.c
Author      :
Date        : 2025-05-11
Description : RSA4096 encryption/decryption implementation.
              This version adds optimized public encryption for e = 65537,
              64-bit Montgomery arithmetic, and CRT-based private decryption.
*****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <stdint.h>

#include "rsa.h"
#include "bignum.h"

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#define RSA_THREAD_LOCAL _Thread_local
#else
#define RSA_THREAD_LOCAL __thread
#endif

#define RSA_FAST_WORDS (RSA_MAX_MODULUS_LEN / 8)

typedef struct {
    int valid;
    uint64_t n[RSA_FAST_WORDS];
    uint64_t rr[RSA_FAST_WORDS];
    uint64_t n0_inv;
} rsa_fast_public_ctx_t;

static rsa_fast_public_ctx_t g_fast_public_ctx;

static int private_block_operation(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_sk_t *sk);
static int public_block_operation_fast_e65537(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_pk_t *pk);

/* Decode big-endian bytes into little-endian 64-bit words. */
static void decode_words64(uint64_t out[RSA_FAST_WORDS], uint8_t *in, uint32_t in_len)
{
    uint32_t i, j;

    for(i=0; i<RSA_FAST_WORDS; i++) {
        uint64_t word = 0;
        for(j=0; j<8; j++) {
            uint32_t pos = in_len - 1 - (i * 8 + j);
            word |= ((uint64_t)in[pos]) << (j * 8);
        }
        out[i] = word;
    }
}

/* Encode little-endian 64-bit words back into big-endian bytes. */
static void encode_words64(uint8_t *out, uint64_t in[RSA_FAST_WORDS])
{
    uint32_t i, j;

    for(i=0; i<RSA_FAST_WORDS; i++) {
        uint64_t word = in[i];
        for(j=0; j<8; j++) {
            out[RSA_MAX_MODULUS_LEN - 1 - (i * 8 + j)] = (uint8_t)(word >> (j * 8));
        }
    }
}

/* Compare two 4096-bit integers stored as 64-bit word arrays. */
static int cmp_words64(uint64_t a[RSA_FAST_WORDS], uint64_t b[RSA_FAST_WORDS])
{
    int i;

    for(i=RSA_FAST_WORDS-1; i>=0; i--) {
        if(a[i] > b[i]) return 1;
        if(a[i] < b[i]) return -1;
    }
    return 0;
}

/* Compute a = a - b for word arrays; caller guarantees a >= b. */
static void sub_words64(uint64_t a[RSA_FAST_WORDS], uint64_t b[RSA_FAST_WORDS])
{
    uint64_t borrow = 0;
    uint32_t i;

    for(i=0; i<RSA_FAST_WORDS; i++) {
        uint64_t bi = b[i] + borrow;
        uint64_t next_borrow = (bi < borrow) || (a[i] < bi);
        a[i] -= bi;
        borrow = next_borrow;
    }
}

/* Compute -n^{-1} mod 2^64 for Montgomery reduction. */
static uint64_t mont_inv64(uint64_t n0)
{
    uint64_t x = 1;
    uint32_t i;

    for(i=0; i<6; i++) {
        x *= (uint64_t)(2 - n0 * x);
    }
    return (uint64_t)(0 - x);
}

/* Compute a = (2 * a) mod n, used to build R^2 mod n. */
static void shift_left_mod64(uint64_t a[RSA_FAST_WORDS], uint64_t n[RSA_FAST_WORDS])
{
    uint64_t carry = 0;
    uint32_t i;

    for(i=0; i<RSA_FAST_WORDS; i++) {
        uint64_t next_carry = a[i] >> 63;
        a[i] = (a[i] << 1) | carry;
        carry = next_carry;
    }

    if(carry || cmp_words64(a, n) >= 0) {
        sub_words64(a, n);
    }
}

/* Compute R^2 mod n for the 4096-bit Montgomery domain. */
static void rr_mont64(uint64_t out[RSA_FAST_WORDS], uint64_t n[RSA_FAST_WORDS])
{
    uint32_t i;

    memset(out, 0, RSA_MAX_MODULUS_LEN);
    out[0] = 1;
    for(i=0; i<2 * RSA_MAX_MODULUS_BITS; i++) {
        shift_left_mod64(out, n);
    }
}

/* 64-bit Montgomery multiplication: out = a * b * R^{-1} mod n. */
static void mont_mul64(uint64_t out[RSA_FAST_WORDS],
                       uint64_t a[RSA_FAST_WORDS],
                       uint64_t b[RSA_FAST_WORDS],
                       uint64_t n[RSA_FAST_WORDS],
                       uint64_t n0_inv)
{
    uint64_t t[RSA_FAST_WORDS + 2] = {0};
    uint32_t i, j;

    for(i=0; i<RSA_FAST_WORDS; i++) {
        __uint128_t carry = 0;
        for(j=0; j<RSA_FAST_WORDS; j++) {
            __uint128_t uv = (__uint128_t)t[j] + (__uint128_t)a[j] * b[i] + carry;
            t[j] = (uint64_t)uv;
            carry = uv >> 64;
        }
        carry += t[RSA_FAST_WORDS];
        t[RSA_FAST_WORDS] = (uint64_t)carry;
        t[RSA_FAST_WORDS + 1] += (uint64_t)(carry >> 64);

        uint64_t m = t[0] * n0_inv;
        carry = 0;
        for(j=0; j<RSA_FAST_WORDS; j++) {
            __uint128_t uv = (__uint128_t)t[j] + (__uint128_t)m * n[j] + carry;
            t[j] = (uint64_t)uv;
            carry = uv >> 64;
        }
        carry += t[RSA_FAST_WORDS];
        t[RSA_FAST_WORDS] = (uint64_t)carry;
        t[RSA_FAST_WORDS + 1] += (uint64_t)(carry >> 64);

        for(j=0; j<=RSA_FAST_WORDS; j++) {
            t[j] = t[j + 1];
        }
        t[RSA_FAST_WORDS + 1] = 0;
    }

    memcpy(out, t, RSA_MAX_MODULUS_LEN);
    if(t[RSA_FAST_WORDS] || cmp_words64(out, n) >= 0) {
        sub_words64(out, n);
    }
}

/*
 * Checks whether the public exponent is 65537.
 * This allows the encryption path to use the optimized e = 65537 routine.
 */
static int pk_exponent_is_65537(rsa_pk_t *pk)
{
    uint32_t i;

    for(i=0; i<RSA_MAX_MODULUS_LEN-3; i++) {
        if(pk->exponent[i] != 0) {
            return 0;
        }
    }

    return pk->exponent[RSA_MAX_MODULUS_LEN-3] == 0x01 &&
           pk->exponent[RSA_MAX_MODULUS_LEN-2] == 0x00 &&
           pk->exponent[RSA_MAX_MODULUS_LEN-1] == 0x01;
}


/*
 * Builds or reuses cached Montgomery parameters for public encryption.
 * The cache avoids recomputing n, R^2 mod n, and n0_inv for every block.
 */
static int get_fast_public_ctx(rsa_pk_t *pk,
                               uint64_t n[RSA_FAST_WORDS],
                               uint64_t rr[RSA_FAST_WORDS],
                               uint64_t *n0_inv)
{
    uint64_t decoded_n[RSA_FAST_WORDS] = {0};

    decode_words64(decoded_n, pk->modulus, RSA_MAX_MODULUS_LEN);
    if((decoded_n[0] & 1) == 0 || ((decoded_n[RSA_FAST_WORDS - 1] >> 63) == 0)) {
        return ERR_WRONG_LEN;
    }

    if(g_fast_public_ctx.valid &&
       cmp_words64(g_fast_public_ctx.n, decoded_n) == 0) {
        memcpy(n, g_fast_public_ctx.n, RSA_MAX_MODULUS_LEN);
        memcpy(rr, g_fast_public_ctx.rr, RSA_MAX_MODULUS_LEN);
        *n0_inv = g_fast_public_ctx.n0_inv;
        return 0;
    }

    memcpy(g_fast_public_ctx.n, decoded_n, RSA_MAX_MODULUS_LEN);
    g_fast_public_ctx.n0_inv = mont_inv64(decoded_n[0]);
    rr_mont64(g_fast_public_ctx.rr, decoded_n);
    g_fast_public_ctx.valid = 1;

    memcpy(n, g_fast_public_ctx.n, RSA_MAX_MODULUS_LEN);
    memcpy(rr, g_fast_public_ctx.rr, RSA_MAX_MODULUS_LEN);
    *n0_inv = g_fast_public_ctx.n0_inv;
    return 0;
}


/*
 * Fast RSA public block operation for e = 65537.
 * Uses 64-bit Montgomery multiplication and computes m^65537 mod n
 * as 16 squarings plus 1 multiplication.
 */
static int public_block_operation_fast_e65537(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_pk_t *pk)
{
    uint64_t m[RSA_FAST_WORDS] = {0};
    uint64_t n[RSA_FAST_WORDS] = {0};
    uint64_t rr[RSA_FAST_WORDS] = {0};
    uint64_t one[RSA_FAST_WORDS] = {0};
    uint64_t base[RSA_FAST_WORDS] = {0};
    uint64_t result[RSA_FAST_WORDS] = {0};
    uint64_t n0_inv;
    uint32_t i;

    if(pk->bits != RSA_MAX_MODULUS_BITS ||
       in_len != RSA_MAX_MODULUS_LEN ||
       !pk_exponent_is_65537(pk)) {
        return ERR_WRONG_LEN;
    }

    decode_words64(m, in, in_len);
    if(get_fast_public_ctx(pk, n, rr, &n0_inv) != 0) {
        return ERR_WRONG_LEN;
    }
    if(cmp_words64(m, n) >= 0) {
        return ERR_WRONG_DATA;
    }

    one[0] = 1;
    mont_mul64(base, m, rr, n, n0_inv);
    memcpy(result, base, RSA_MAX_MODULUS_LEN);

    for(i=0; i<16; i++) {
        mont_mul64(result, result, result, n, n0_inv);
    }
    mont_mul64(result, result, base, n, n0_inv);
    mont_mul64(result, result, one, n, n0_inv);

    *out_len = RSA_MAX_MODULUS_LEN;
    encode_words64(out, result);

    memset(m, 0, sizeof(m));
    memset(n, 0, sizeof(n));
    memset(rr, 0, sizeof(rr));
    memset(base, 0, sizeof(base));
    memset(result, 0, sizeof(result));

    return 0;
}

int rsa_private_encrypt_any_len(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_sk_t *sk){
	int status=0;
	uint32_t len=0;
	uint8_t *tmp_o=out;
	*out_len=0;
	for(uint32_t i=0;i<in_len && status==0;i+=(RSA_MAX_MODULUS_LEN-11)){
		if((in_len-i)>(RSA_MAX_MODULUS_LEN-11)){
			status=rsa_private_encrypt(tmp_o,&len,in+i,RSA_MAX_MODULUS_LEN-11,sk);
		}
		else{
			status=rsa_private_encrypt(tmp_o,&len,in+i,in_len-i,sk);
			*out_len+=len;
			break;
		}
		tmp_o=tmp_o+len;
		*out_len+=len;
	}
	tmp_o=NULL;
	return status;
}

/* Encrypt arbitrary-length input by splitting it into PKCS#1-sized blocks. */
int rsa_public_encrypt_any_len(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_pk_t *pk){
    uint32_t modulus_len = (pk->bits + 7) / 8;
    uint32_t max_plain_len = modulus_len - 11;
    uint8_t *tmp_o = out;
    uint32_t len = 0;
    int status = 0;
    uint32_t off;

    *out_len = 0;
    for(off=0; off<in_len && status==0; off+=max_plain_len) {
        uint32_t chunk_len = in_len - off;
        if(chunk_len > max_plain_len) {
            chunk_len = max_plain_len;
        }
        status = rsa_public_encrypt(tmp_o, &len, in + off, chunk_len, pk);
        tmp_o += len;
        *out_len += len;
    }
    return status;
}


/* Decrypt arbitrary-length ciphertext by processing one RSA block at a time. */
int rsa_private_decrypt_any_len(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_sk_t *sk){
    uint32_t modulus_len = (sk->bits + 7) / 8;
    uint8_t *tmp_o = out;
    uint32_t len = 0;
    int status = 0;
    uint32_t off;

    *out_len = 0;
    for(off=0; off<in_len && status==0; off+=modulus_len) {
        uint32_t chunk_len = in_len - off;
        if(chunk_len > modulus_len) {
            chunk_len = modulus_len;
        }
        status = rsa_private_decrypt(tmp_o, &len, in + off, chunk_len, sk);
        tmp_o += len;
        *out_len += len;
    }
    return status;
}

int rsa_private_encrypt(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_sk_t *sk)
{
    int status;
    uint8_t pkcs_block[RSA_MAX_MODULUS_LEN];
    uint32_t i, modulus_len;

    modulus_len = (sk->bits + 7) / 8;
    if(in_len + 11 > modulus_len)
        return ERR_WRONG_LEN;

    pkcs_block[0] = 0;
    pkcs_block[1] = 1;
    for(i=2; i<modulus_len-in_len-1; i++) {
        pkcs_block[i] = 0xFF;
    }

    pkcs_block[i++] = 0;

    memcpy((uint8_t *)&pkcs_block[i], (uint8_t *)in, in_len);

    status = private_block_operation(out, out_len, pkcs_block, modulus_len, sk);

    // Clear potentially sensitive information
    memset((uint8_t *)pkcs_block, 0, sizeof(pkcs_block));

    return status;
}

int rsa_private_decrypt(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_sk_t *sk)
{
    int status;
    uint8_t pkcs_block[RSA_MAX_MODULUS_LEN];
    uint32_t i, modulus_len, pkcs_block_len;

    modulus_len = (sk->bits + 7) / 8;
    if(in_len > modulus_len)
        return ERR_WRONG_LEN;

    status = private_block_operation(pkcs_block, &pkcs_block_len, in, in_len, sk);
    if(status != 0)
        return status;

    if(pkcs_block_len != modulus_len)
        return ERR_WRONG_LEN;

    if((pkcs_block[0] != 0) || (pkcs_block[1] != 2))
        return ERR_WRONG_DATA;

    for(i=2; i<modulus_len-1; i++) {
        if(pkcs_block[i] == 0)  break;
    }

    i++;
    if(i >= modulus_len)
        return ERR_WRONG_DATA;
    *out_len = modulus_len - i;
    if(*out_len + 11 > modulus_len)
        return ERR_WRONG_DATA;
    memcpy((uint8_t *)out, (uint8_t *)&pkcs_block[i], *out_len);
    // Clear potentially sensitive information
    memset((uint8_t *)pkcs_block, 0, sizeof(pkcs_block));

    return status;
}

static int private_block_operation(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_sk_t *sk)
{
    uint32_t ndigits, pdigits, qdigits, dpdigits, dqdigits, prime_digits;
    bn_t c[BN_MAX_DIGITS] = {0}, n[BN_MAX_DIGITS] = {0};
    bn_t p[BN_MAX_DIGITS] = {0}, q[BN_MAX_DIGITS] = {0};
    bn_t dp[BN_MAX_DIGITS] = {0}, dq[BN_MAX_DIGITS] = {0}, qinv[BN_MAX_DIGITS] = {0};
    bn_t cp[BN_MAX_DIGITS] = {0}, cq[BN_MAX_DIGITS] = {0};
    bn_t m1[BN_MAX_DIGITS] = {0}, m2[BN_MAX_DIGITS] = {0}, m2_mod_p[BN_MAX_DIGITS] = {0};
    bn_t diff[BN_MAX_DIGITS] = {0}, h[BN_MAX_DIGITS] = {0};
    bn_t qh[2 * BN_MAX_DIGITS] = {0}, m2_wide[2 * BN_MAX_DIGITS] = {0};
    bn_t m[2 * BN_MAX_DIGITS] = {0};

    bn_decode(c, BN_MAX_DIGITS, in, in_len);
    bn_decode(n, BN_MAX_DIGITS, sk->modulus, RSA_MAX_MODULUS_LEN);
    bn_decode(p, BN_MAX_DIGITS, sk->prime1, RSA_MAX_PRIME_LEN);
    bn_decode(q, BN_MAX_DIGITS, sk->prime2, RSA_MAX_PRIME_LEN);
    bn_decode(dp, BN_MAX_DIGITS, sk->prime_exponent1, RSA_MAX_PRIME_LEN);
    bn_decode(dq, BN_MAX_DIGITS, sk->prime_exponent2, RSA_MAX_PRIME_LEN);
    bn_decode(qinv, BN_MAX_DIGITS, sk->coefficient, RSA_MAX_PRIME_LEN);

    ndigits = bn_digits(n, BN_MAX_DIGITS);
    pdigits = bn_digits(p, BN_MAX_DIGITS);
    qdigits = bn_digits(q, BN_MAX_DIGITS);
    dpdigits = bn_digits(dp, BN_MAX_DIGITS);
    dqdigits = bn_digits(dq, BN_MAX_DIGITS);
    prime_digits = (pdigits > qdigits) ? pdigits : qdigits;

    if(bn_cmp(c, n, ndigits) >= 0)
        return ERR_WRONG_DATA;

    if(pdigits == 0 || qdigits == 0 || dpdigits == 0 || dqdigits == 0)
        return ERR_WRONG_DATA;

    // CRT: m1 = c^dp mod p, m2 = c^dq mod q, h = qInv * (m1 - m2) mod p.
    bn_mod(cp, c, ndigits, p, pdigits);
    bn_mod(cq, c, ndigits, q, qdigits);
    bn_mod_exp_mont(m1, cp, dp, dpdigits, p, pdigits);
    bn_mod_exp_mont(m2, cq, dq, dqdigits, q, qdigits);

    bn_mod(m2_mod_p, m2, qdigits, p, pdigits);
    if(bn_sub(diff, m1, m2_mod_p, pdigits)) {
        bn_add(diff, diff, p, pdigits);
    }
    bn_mod_mul_mont(h, qinv, diff, p, pdigits);

    bn_mul(qh, q, h, prime_digits);
    bn_assign(m2_wide, m2, qdigits);
    bn_add(m, qh, m2_wide, 2 * prime_digits);

    *out_len = (sk->bits + 7) / 8;
    bn_encode(out, *out_len, m, ndigits);

    // Clear potentially sensitive information
    memset((uint8_t *)c, 0, sizeof(c));
    memset((uint8_t *)n, 0, sizeof(n));
    memset((uint8_t *)p, 0, sizeof(p));
    memset((uint8_t *)q, 0, sizeof(q));
    memset((uint8_t *)dp, 0, sizeof(dp));
    memset((uint8_t *)dq, 0, sizeof(dq));
    memset((uint8_t *)qinv, 0, sizeof(qinv));
    memset((uint8_t *)cp, 0, sizeof(cp));
    memset((uint8_t *)cq, 0, sizeof(cq));
    memset((uint8_t *)m1, 0, sizeof(m1));
    memset((uint8_t *)m2, 0, sizeof(m2));
    memset((uint8_t *)m2_mod_p, 0, sizeof(m2_mod_p));
    memset((uint8_t *)diff, 0, sizeof(diff));
    memset((uint8_t *)h, 0, sizeof(h));
    memset((uint8_t *)qh, 0, sizeof(qh));
    memset((uint8_t *)m2_wide, 0, sizeof(m2_wide));
    memset((uint8_t *)m, 0, sizeof(m));

    return 0;
}

// Public encryption
static int public_block_operation(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_pk_t *pk);

void generate_rand(uint8_t *block, uint32_t block_len)
{
    uint32_t i;
    static RSA_THREAD_LOCAL uint32_t seed = 0;
    if(seed == 0) {
        seed = (uint32_t)time(NULL) ^ (uint32_t)(uintptr_t)&seed;
        if(seed == 0) {
            seed = 1;
        }
    }
    for(i=0; i<block_len; i++) {
        do {
            seed = seed * 1664525u + 1013904223u;
            block[i] = (uint8_t)(seed >> 24);
        } while(block[i] == 0);
    }
}

int rsa_public_encrypt(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_pk_t *pk)
{
    int status;
    uint8_t byte, pkcs_block[RSA_MAX_MODULUS_LEN];
    uint32_t i, modulus_len;

    modulus_len = (pk->bits + 7) / 8;
    if(in_len + 11 > modulus_len) {//padding len
        return ERR_WRONG_LEN;
    }

    pkcs_block[0] = 0;
    pkcs_block[1] = 2;
    for(i=2; i<modulus_len-in_len-1; i++) {
        do {
            generate_rand(&byte, 1);
        } while(byte == 0);
        pkcs_block[i] = byte;
    }
    pkcs_block[i++] = 0;

    memcpy((uint8_t *)&pkcs_block[i], (uint8_t *)in, in_len);
    status = public_block_operation(out, out_len, pkcs_block, modulus_len, pk);
    // Clear potentially sensitive information
    byte = 0;
    memset((uint8_t *)pkcs_block, 0, sizeof(pkcs_block));

    return status;
}



static int public_block_operation(uint8_t *out, uint32_t *out_len, uint8_t *in, uint32_t in_len, rsa_pk_t *pk)
{
    int fast_status;
    uint32_t edigits, ndigits;
    bn_t c[BN_MAX_DIGITS], e[BN_MAX_DIGITS], m[BN_MAX_DIGITS], n[BN_MAX_DIGITS];

    fast_status = public_block_operation_fast_e65537(out, out_len, in, in_len, pk);
    if(fast_status == 0 || fast_status == ERR_WRONG_DATA) {
        return fast_status;
    }

    bn_decode(m, BN_MAX_DIGITS, in, in_len);
    bn_decode(n, BN_MAX_DIGITS, pk->modulus, RSA_MAX_MODULUS_LEN);
    bn_decode(e, BN_MAX_DIGITS, pk->exponent, RSA_MAX_MODULUS_LEN);

    ndigits = bn_digits(n, BN_MAX_DIGITS);
    edigits = bn_digits(e, BN_MAX_DIGITS);

    if(bn_cmp(m, n, ndigits) >= 0) {
        return ERR_WRONG_DATA;
    }

    if(edigits == 1 && e[0] == 0x00010001) {
        bn_mod_exp_mont_e65537(c, m, n, ndigits);
    } else {
        bn_mod_exp_mont(c, m, e, edigits, n, ndigits);
    }

    *out_len = (pk->bits + 7) / 8;
    bn_encode(out, *out_len, c, ndigits);

    // Clear potentially sensitive information
    memset((uint8_t *)c, 0, sizeof(c));
    memset((uint8_t *)m, 0, sizeof(m));

    return 0;
}
