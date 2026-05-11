/*****************************************************************************
Filename    : bignum.c
Author      :
Date        : 2025-5-11
Description : core file for big number arithmetic operations.
*****************************************************************************/
#include <string.h>
#include "bignum.h"

#define BN_MONT_CACHE_SLOTS 4
#define BN_MONT_WINDOW_BITS 6
#define BN_MONT_TABLE_SIZE (1u << (BN_MONT_WINDOW_BITS - 1))

typedef struct {
    int valid;
    uint32_t digits;
    bn_t modulus[BN_MAX_DIGITS];
    bn_t rr[BN_MAX_DIGITS];
    bn_t one_mont[BN_MAX_DIGITS];
    bn_t n0_inv;
} bn_mont_cache_t;

static bn_mont_cache_t bn_mont_cache[BN_MONT_CACHE_SLOTS];
static uint32_t bn_mont_cache_next = 0;

static bn_t bn_sub_digit_mul(bn_t *a, bn_t *b, bn_t c, bn_t *d, uint32_t digits);
static bn_t bn_add_digit_mul(bn_t *a, bn_t *b, bn_t c, bn_t *d, uint32_t digits);
static uint32_t bn_digit_bits(bn_t a);
static bn_t bn_mont_inv32(bn_t n0);
static uint32_t bn_get_window(bn_t *a, uint32_t bit, uint32_t width);
static void bn_mont_shift_mod(bn_t *a, bn_t *n, uint32_t digits, uint32_t shift_digits);
static void bn_get_mont_context(bn_t *rr, bn_t *one_mont, bn_t *n0_inv, bn_t *n, uint32_t digits);
static void bn_mont_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *n, uint32_t digits, bn_t n0_inv);

void bn_decode(bn_t *bn, uint32_t digits, uint8_t *hexarr, uint32_t size)
{
    bn_t t;
    int j;
    uint32_t i, u;
    for(i=0,j=size-1; i<digits && j>=0; i++) {
        t = 0;
        for(u=0; j>=0 && u<BN_DIGIT_BITS; j--, u+=8) {
            t |= ((bn_t)hexarr[j]) << u;
        }
        bn[i] = t;
    }

    for(; i<digits; i++) {
        bn[i] = 0;
    }
}

void bn_encode(uint8_t *hexarr, uint32_t size, bn_t *bn, uint32_t digits)
{
    bn_t t;
    int j;
    uint32_t i, u;

    for(i=0,j=size-1; i<digits && j>=0; i++) {
        t = bn[i];
        for(u=0; j>=0 && u<BN_DIGIT_BITS; j--, u+=8) {
            hexarr[j] = (uint8_t)(t >> u);
        }
    }

    for(; j>=0; j--) {
        hexarr[j] = 0;
    }
}

void bn_assign(bn_t *a, bn_t *b, uint32_t digits)
{
    uint32_t i;
    for(i=0; i<digits; i++) {
        a[i] = b[i];
    }
}

void bn_assign_zero(bn_t *a, uint32_t digits)
{
    uint32_t i;
    for(i=0; i<digits; i++) {
        a[i] = 0;
    }
}

bn_t bn_add(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    bn_t ai, carry;
    uint32_t i;

    carry = 0;
    for(i=0; i<digits; i++) {
        if((ai = b[i] + carry) < carry) {
            ai = c[i];
        } else if((ai += c[i]) < c[i]) {
            carry = 1;
        } else {
            carry = 0;
        }
        a[i] = ai;
    }

    return carry;
}

bn_t bn_sub(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    bn_t ai, borrow;
    uint32_t i;

    borrow = 0;
    for(i=0; i<digits; i++) {
        if((ai = b[i] - borrow) > (BN_MAX_DIGIT - borrow)) {
            ai = BN_MAX_DIGIT - c[i];
        } else if((ai -= c[i]) > (BN_MAX_DIGIT - c[i])) {
            borrow = 1;
        } else {
            borrow = 0;
        }
        a[i] = ai;
    }

    return borrow;
}

void bn_mul(bn_t *a, bn_t *b, bn_t *c, uint32_t digits)
{
    bn_t t[2*BN_MAX_DIGITS];
    uint32_t bdigits, cdigits, i;

    bn_assign_zero(t, 2*digits);
    bdigits = bn_digits(b, digits);
    cdigits = bn_digits(c, digits);

    for(i=0; i<bdigits; i++) {
        t[i+cdigits] += bn_add_digit_mul(&t[i], &t[i], b[i], c, cdigits);
    }

    bn_assign(a, t, 2*digits);

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
}

void bn_div(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *d, uint32_t ddigits)
{
    dbn_t tmp;
    bn_t ai, t, cc[2*BN_MAX_DIGITS+1], dd[BN_MAX_DIGITS];
    int i;
    uint32_t dddigits, shift;

    dddigits = bn_digits(d, ddigits);
    if(dddigits == 0)
        return;

    shift = BN_DIGIT_BITS - bn_digit_bits(d[dddigits-1]);
    bn_assign_zero(cc, dddigits);
    cc[cdigits] = bn_shift_l(cc, c, shift, cdigits);
    bn_shift_l(dd, d, shift, dddigits);
    t = dd[dddigits-1];

    bn_assign_zero(a, cdigits);
    i = cdigits - dddigits;
    for(; i>=0; i--) {
        if(t == BN_MAX_DIGIT) {
            ai = cc[i+dddigits];
        } else {
            tmp = cc[i+dddigits-1];
            tmp += (dbn_t)cc[i+dddigits] << BN_DIGIT_BITS;
            ai = tmp / (t + 1);
        }

        cc[i+dddigits] -= bn_sub_digit_mul(&cc[i], &cc[i], ai, dd, dddigits);
        // printf("cc[%d]: %08X\n", i, cc[i+dddigits]);
        while(cc[i+dddigits] || (bn_cmp(&cc[i], dd, dddigits) >= 0)) {
            ai++;
            cc[i+dddigits] -= bn_sub(&cc[i], &cc[i], dd, dddigits);
        }
        a[i] = ai;
    }

    bn_assign_zero(b, ddigits);
    bn_shift_r(b, cc, shift, dddigits);

    // Clear potentially sensitive information
    memset((uint8_t *)cc, 0, sizeof(cc));
    memset((uint8_t *)dd, 0, sizeof(dd));
}

bn_t bn_shift_l(bn_t *a, bn_t *b, uint32_t c, uint32_t digits)
{
    bn_t bi, carry;
    uint32_t i, t;

    if(c >= BN_DIGIT_BITS)
        return 0;

    t = BN_DIGIT_BITS - c;
    carry = 0;
    for(i=0; i<digits; i++) {
        bi = b[i];
        a[i] = (bi << c) | carry;
        carry = c ? (bi >> t) : 0;
    }

    return carry;
}

bn_t bn_shift_r(bn_t *a, bn_t *b, uint32_t c, uint32_t digits)
{
    bn_t bi, carry;
    int i;
    uint32_t t;

    if(c >= BN_DIGIT_BITS)
        return 0;

    t = BN_DIGIT_BITS - c;
    carry = 0;
    i = digits - 1;
    for(; i>=0; i--) {
        bi = b[i];
        a[i] = (bi >> c) | carry;
        carry = c ? (bi << t) : 0;
    }

    return carry;
}

void bn_mod(bn_t *a, bn_t *b, uint32_t bdigits, bn_t *c, uint32_t cdigits)
{
    bn_t t[2*BN_MAX_DIGITS] = {0};

    bn_div(t, a, b, bdigits, c, cdigits);

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
}

void bn_mod_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *d, uint32_t digits)
{
    bn_t t[2*BN_MAX_DIGITS];

    bn_mul(t, b, c, digits);
    bn_mod(a, t, 2*digits, d, digits);

    // Clear potentially sensitive information
    memset((uint8_t *)t, 0, sizeof(t));
}

/*
 * Montgomery modular multiplication wrapper.
 * Computes a = b * c mod d.
 * For odd modulus d, operands are converted into Montgomery form, multiplied,
 * and converted back. For even modulus d, it falls back to bn_mod_mul().
 */
void bn_mod_mul_mont(bn_t *a, bn_t *b, bn_t *c, bn_t *d, uint32_t digits)
{
    bn_t bb[BN_MAX_DIGITS] = {0}, cc[BN_MAX_DIGITS] = {0};
    bn_t bm[BN_MAX_DIGITS] = {0}, cm[BN_MAX_DIGITS] = {0}, rm[BN_MAX_DIGITS] = {0};
    bn_t rr[BN_MAX_DIGITS] = {0}, mont_tmp[BN_MAX_DIGITS] = {0}, one[BN_MAX_DIGITS] = {0};
    bn_t n0_inv;

    if(digits == 0)
        return;
    if((d[0] & 1) == 0) {
        bn_mod_mul(a, b, c, d, digits);
        return;
    }

    bn_assign(bb, b, digits);
    bn_assign(cc, c, digits);
    if(bn_cmp(bb, d, digits) >= 0) {
        bn_mod(bb, bb, digits, d, digits);
    }
    if(bn_cmp(cc, d, digits) >= 0) {
        bn_mod(cc, cc, digits, d, digits);
    }

    one[0] = 1;
    bn_get_mont_context(rr, mont_tmp, &n0_inv, d, digits);
    bn_mont_mul(bm, bb, rr, d, digits, n0_inv);
    bn_mont_mul(cm, cc, rr, d, digits, n0_inv);
    bn_mont_mul(rm, bm, cm, d, digits, n0_inv);
    bn_mont_mul(a, rm, one, d, digits, n0_inv);

    memset((uint8_t *)bb, 0, sizeof(bb));
    memset((uint8_t *)cc, 0, sizeof(cc));
    memset((uint8_t *)bm, 0, sizeof(bm));
    memset((uint8_t *)cm, 0, sizeof(cm));
    memset((uint8_t *)rm, 0, sizeof(rm));
    memset((uint8_t *)rr, 0, sizeof(rr));
    memset((uint8_t *)mont_tmp, 0, sizeof(mont_tmp));
    memset((uint8_t *)one, 0, sizeof(one));
}


/*
 * Montgomery modular exponentiation with sliding-window method.
 * Computes a = b^c mod d.
 * Precomputes odd powers of the base in Montgomery form to reduce the number
 * of modular multiplications during RSA exponentiation.
 */
void bn_mod_exp_mont(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *d, uint32_t ddigits)
{
    bn_t base[BN_MAX_DIGITS] = {0}, base2[BN_MAX_DIGITS] = {0};
    bn_t result[BN_MAX_DIGITS] = {0}, one_mont[BN_MAX_DIGITS] = {0};
    bn_t rr[BN_MAX_DIGITS] = {0}, one[BN_MAX_DIGITS] = {0};
    bn_t table[BN_MONT_TABLE_SIZE][BN_MAX_DIGITS] = {{0}};
    bn_t n0_inv;
    uint32_t i, total_bits;
    int bit;

    if(ddigits == 0)
        return;
    if((d[0] & 1) == 0) {
        bn_mod_exp(a, b, c, cdigits, d, ddigits);
        return;
    }

    one[0] = 1;
    bn_get_mont_context(rr, one_mont, &n0_inv, d, ddigits);
    bn_mont_mul(base, b, rr, d, ddigits, n0_inv);

    cdigits = bn_digits(c, cdigits);
    if(cdigits == 0) {
        bn_assign(a, one, ddigits);
        return;
    }

    bn_assign(table[0], base, ddigits);
    bn_mont_mul(base2, base, base, d, ddigits, n0_inv);
    for(i=1; i<BN_MONT_TABLE_SIZE; i++) {
        bn_mont_mul(table[i], table[i - 1], base2, d, ddigits, n0_inv);
    }

    bn_assign(result, one_mont, ddigits);
    total_bits = (cdigits - 1) * BN_DIGIT_BITS + bn_digit_bits(c[cdigits - 1]);
    bit = (int)total_bits - 1;
    while(bit >= 0) {
        uint32_t width, start, w;

        if(bn_get_window(c, (uint32_t)bit, 1) == 0) {
            bn_mont_mul(result, result, result, d, ddigits, n0_inv);
            bit--;
            continue;
        }

        width = BN_MONT_WINDOW_BITS;
        if(width > (uint32_t)bit + 1) {
            width = (uint32_t)bit + 1;
        }
        start = (uint32_t)bit + 1 - width;
        w = bn_get_window(c, start, width);
        while((w & 1u) == 0) {
            w >>= 1;
            width--;
            start++;
        }

        for(i=0; i<width; i++) {
            bn_mont_mul(result, result, result, d, ddigits, n0_inv);
        }
        bn_mont_mul(result, result, table[w >> 1], d, ddigits, n0_inv);
        bit = (int)start - 1;
    }

    bn_mont_mul(a, result, one, d, ddigits, n0_inv);

    memset((uint8_t *)base, 0, sizeof(base));
    memset((uint8_t *)base2, 0, sizeof(base2));
    memset((uint8_t *)result, 0, sizeof(result));
    memset((uint8_t *)one_mont, 0, sizeof(one_mont));
    memset((uint8_t *)rr, 0, sizeof(rr));
    memset((uint8_t *)one, 0, sizeof(one));
    memset((uint8_t *)table, 0, sizeof(table));
}


/*
 * Fast modular exponentiation for RSA public exponent e = 65537.
 * Since 65537 = 2^16 + 1, encryption only needs 16 modular squarings
 * and 1 modular multiplication in Montgomery form.
 */
void bn_mod_exp_mont_e65537(bn_t *a, bn_t *b, bn_t *d, uint32_t ddigits)
{
    bn_t base[BN_MAX_DIGITS] = {0}, result[BN_MAX_DIGITS] = {0};
    bn_t rr[BN_MAX_DIGITS] = {0}, one_mont[BN_MAX_DIGITS] = {0};
    bn_t one[BN_MAX_DIGITS] = {0};
    bn_t n0_inv;
    uint32_t i;

    if(ddigits == 0)
        return;
    if((d[0] & 1) == 0) {
        bn_assign(result, b, ddigits);
        for(i=0; i<16; i++) {
            bn_mod_mul(result, result, result, d, ddigits);
        }
        bn_mod_mul(a, result, b, d, ddigits);
        return;
    }

    one[0] = 1;
    bn_get_mont_context(rr, one_mont, &n0_inv, d, ddigits);
    bn_mont_mul(base, b, rr, d, ddigits, n0_inv);
    bn_assign(result, base, ddigits);
    for(i=0; i<16; i++) {
        bn_mont_mul(result, result, result, d, ddigits, n0_inv);
    }
    bn_mont_mul(result, result, base, d, ddigits, n0_inv);
    bn_mont_mul(a, result, one, d, ddigits, n0_inv);

    memset((uint8_t *)base, 0, sizeof(base));
    memset((uint8_t *)result, 0, sizeof(result));
    memset((uint8_t *)rr, 0, sizeof(rr));
    memset((uint8_t *)one_mont, 0, sizeof(one_mont));
    memset((uint8_t *)one, 0, sizeof(one));
}



void bn_mod_exp(bn_t *a, bn_t *b, bn_t *c, uint32_t cdigits, bn_t *d, uint32_t ddigits)
{
    bn_t bpower[3][BN_MAX_DIGITS], ci, t[BN_MAX_DIGITS];
    int i;
    uint32_t ci_bits, j, s;

    bn_assign(bpower[0], b, ddigits);
    bn_mod_mul(bpower[1], bpower[0], b, d, ddigits);
    bn_mod_mul(bpower[2], bpower[1], b, d, ddigits);

    BN_ASSIGN_DIGIT(t, 1, ddigits);

    cdigits = bn_digits(c, cdigits);
    i = cdigits - 1;
    for(; i>=0; i--) {
        ci = c[i];
        ci_bits = BN_DIGIT_BITS;

        if(i == (int)(cdigits - 1)) {
            while(!DIGIT_2MSB(ci)) {
                ci <<= 2;
                ci_bits -= 2;
            }
        }

        for(j=0; j<ci_bits; j+=2) {
            bn_mod_mul(t, t, t, d, ddigits);
            bn_mod_mul(t, t, t, d, ddigits);
            if((s = DIGIT_2MSB(ci)) != 0) {
                bn_mod_mul(t, t, bpower[s-1], d, ddigits);
            }
            ci <<= 2;
        }
    }

    bn_assign(a, t, ddigits);

    // Clear potentially sensitive information
    memset((uint8_t *)bpower, 0, sizeof(bpower));
    memset((uint8_t *)t, 0, sizeof(t));
}

int bn_cmp(bn_t *a, bn_t *b, uint32_t digits)
{
    int i;
    for(i=digits-1; i>=0; i--) {
        if(a[i] > b[i])     return 1;
        if(a[i] < b[i])     return -1;
    }

    return 0;
}

uint32_t bn_digits(bn_t *a, uint32_t digits)
{
    int i;
    for(i=digits-1; i>=0; i--) {
        if(a[i])    break;
    }

    return (i + 1);
}

static bn_t bn_add_digit_mul(bn_t *a, bn_t *b, bn_t c, bn_t *d, uint32_t digits)
{
    dbn_t result;
    bn_t carry, rh, rl;
    uint32_t i;

    if(c == 0)
        return 0;

    carry = 0;
    for(i=0; i<digits; i++) {
        result = (dbn_t)c * d[i];
        rl = result & BN_MAX_DIGIT;
        rh = (result >> BN_DIGIT_BITS) & BN_MAX_DIGIT;
        if((a[i] = b[i] + carry) < carry) {
            carry = 1;
        } else {
            carry = 0;
        }
        if((a[i] += rl) < rl) {
            carry++;
        }
        carry += rh;
    }

    return carry;
}

static bn_t bn_sub_digit_mul(bn_t *a, bn_t *b, bn_t c, bn_t *d, uint32_t digits)
{
    dbn_t result;
    bn_t borrow, rh, rl;
    uint32_t i;

    if(c == 0)
        return 0;

    borrow = 0;
    for(i=0; i<digits; i++) {
        result = (dbn_t)c * d[i];
        rl = result & BN_MAX_DIGIT;
        rh = (result >> BN_DIGIT_BITS) & BN_MAX_DIGIT;
        if((a[i] = b[i] - borrow) > (BN_MAX_DIGIT - borrow)) {
            borrow = 1;
        } else {
            borrow = 0;
        }
        if((a[i] -= rl) > (BN_MAX_DIGIT - rl)) {
            borrow++;
        }
        borrow += rh;
    }

    return borrow;
}

static uint32_t bn_digit_bits(bn_t a)
{
    uint32_t i;
    for(i=0; i<BN_DIGIT_BITS; i++) {
        if(a == 0)  break;
        a >>= 1;
    }

    return i;
}

/* Compute -n^{-1} mod 2^32 for Montgomery reduction. */
static bn_t bn_mont_inv32(bn_t n0)
{
    bn_t x = 1;
    uint32_t i;

    for(i=0; i<5; i++) {
        x *= (bn_t)(2 - n0 * x);
    }

    return (bn_t)(0 - x);
}

/* Read a small bit window from a multi-precision exponent. */
static uint32_t bn_get_window(bn_t *a, uint32_t bit, uint32_t width)
{
    uint32_t digit = bit / BN_DIGIT_BITS;
    uint32_t offset = bit % BN_DIGIT_BITS;
    uint64_t value = a[digit] >> offset;
    uint32_t mask = (1u << width) - 1u;

    if(offset + width > BN_DIGIT_BITS && digit + 1 < BN_MAX_DIGITS) {
        value |= (uint64_t)a[digit + 1] << (BN_DIGIT_BITS - offset);
    }

    return (uint32_t)value & mask;
}

/* Compute R or R^2 modulo n for Montgomery representation. */
static void bn_mont_shift_mod(bn_t *a, bn_t *n, uint32_t digits, uint32_t shift_digits)
{
    bn_t shifted[2 * BN_MAX_DIGITS + 1] = {0};

    shifted[shift_digits] = 1;
    bn_mod(a, shifted, shift_digits + 1, n, digits);

    memset((uint8_t *)shifted, 0, sizeof(shifted));
}

/*
 * Build or reuse Montgomery constants for one modulus.
 * Cached values avoid recomputing R mod n, R^2 mod n, and n0_inv repeatedly.
 */
static void bn_get_mont_context(bn_t *rr, bn_t *one_mont, bn_t *n0_inv, bn_t *n, uint32_t digits)
{
    uint32_t i, slot;

    for(i=0; i<BN_MONT_CACHE_SLOTS; i++) {
        if(bn_mont_cache[i].valid &&
           bn_mont_cache[i].digits == digits &&
           bn_cmp(bn_mont_cache[i].modulus, n, digits) == 0) {
            bn_assign(rr, bn_mont_cache[i].rr, digits);
            bn_assign(one_mont, bn_mont_cache[i].one_mont, digits);
            *n0_inv = bn_mont_cache[i].n0_inv;
            return;
        }
    }

    slot = bn_mont_cache_next;
    bn_mont_cache_next = (bn_mont_cache_next + 1) % BN_MONT_CACHE_SLOTS;
    bn_mont_cache[slot].valid = 1;
    bn_mont_cache[slot].digits = digits;
    bn_assign_zero(bn_mont_cache[slot].modulus, BN_MAX_DIGITS);
    bn_assign_zero(bn_mont_cache[slot].rr, BN_MAX_DIGITS);
    bn_assign_zero(bn_mont_cache[slot].one_mont, BN_MAX_DIGITS);
    bn_assign(bn_mont_cache[slot].modulus, n, digits);
    bn_mont_cache[slot].n0_inv = bn_mont_inv32(n[0]);
    bn_mont_shift_mod(bn_mont_cache[slot].one_mont, n, digits, digits);
    bn_mont_shift_mod(bn_mont_cache[slot].rr, n, digits, 2 * digits);

    bn_assign(rr, bn_mont_cache[slot].rr, digits);
    bn_assign(one_mont, bn_mont_cache[slot].one_mont, digits);
    *n0_inv = bn_mont_cache[slot].n0_inv;
}


/*
 * Core Montgomery reduction multiplication.
 * Inputs b and c are Montgomery-domain values.
 * Computes a = b * c * R^{-1} mod n without division by n.
 */
static void bn_mont_mul(bn_t *a, bn_t *b, bn_t *c, bn_t *n, uint32_t digits, bn_t n0_inv)
{
    dbn_t t[BN_MAX_DIGITS + 2] = {0};
    dbn_t uv, carry;
    bn_t m;
    uint32_t i, j;

    for(i=0; i<digits; i++) {
        carry = 0;
        for(j=0; j<digits; j++) {
            uv = t[j] + (dbn_t)b[j] * c[i] + carry;
            t[j] = (bn_t)uv;
            carry = uv >> BN_DIGIT_BITS;
        }
        uv = t[digits] + carry;
        t[digits] = (bn_t)uv;
        t[digits + 1] += uv >> BN_DIGIT_BITS;

        m = (bn_t)t[0] * n0_inv;
        carry = 0;
        uv = t[0] + (dbn_t)m * n[0];
        carry = uv >> BN_DIGIT_BITS;
        for(j=1; j<digits; j++) {
            uv = t[j] + (dbn_t)m * n[j] + carry;
            t[j - 1] = (bn_t)uv;
            carry = uv >> BN_DIGIT_BITS;
        }
        uv = t[digits] + carry;
        t[digits - 1] = (bn_t)uv;
        t[digits] = t[digits + 1] + (uv >> BN_DIGIT_BITS);
        t[digits + 1] = 0;
    }

    for(i=0; i<digits; i++) {
        a[i] = (bn_t)t[i];
    }
    if(t[digits] || bn_cmp(a, n, digits) >= 0) {
        bn_sub(a, a, n, digits);
    }

    memset((uint8_t *)t, 0, sizeof(t));
}
