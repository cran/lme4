#include "bitfield.h"
#include <limits.h>

/**
 * Actual definition of struct nlme_bitfield.
 * len is the number of bits.
 * words contains the bits.
 * 
 */
struct nlme_bitfield {
    nlme_word len;
    nlme_word *words;
};

#define WORD_BITS (CHAR_BIT*sizeof(int))

#ifndef min
#define min(x, y) (((x)<(y))?(x):(y))
#endif

/* #define NLME_TEST_BITFIELD */

#ifndef NLME_TEST_BITFIELD
#include "R.h"
#else
#include <stdio.h>
#include <stdlib.h>

void
error(const char *mes)
{
    fprintf(stderr, mes);
    fprintf(stderr, "\n");
    exit(0);
}

#define Calloc(n, t) ((t *) calloc((n)*sizeof(t)))
#define Free(x) free(x)
#define REprintf printf
#endif

void
nlme_bitfield_print(nlme_bitfield *x, nlme_word start,
                    nlme_word end)
{
    nlme_word i;

    if (x == NULL)
        return;

    if (x->len <= start || x->len < end)
        error("Invalid index to bitfield");

    for (i = (start/WORD_BITS)*WORD_BITS; i < start; i++) {
        REprintf(" ");
    }
    for (i = start; i < end; i++) {
        nlme_word word = x->words[i/WORD_BITS];
        REprintf((word & (1 << (i%WORD_BITS)))?"1":"0");
        if ((i+1) % WORD_BITS == 0)
            REprintf("\n");
    }
    if ((i+1) % WORD_BITS != 0)
        REprintf("\n");
}

#ifdef NLME_TEST_BITFIELD

int
main(int argc, char **argv)
{
    nlme_word len = (argc == 1)?(1 << 16):atoi(argv[1]);
    nlme_range range;
    nlme_bitfield *bf = nlme_bitfield_alloc(400);
    int i;

    nlme_bitfield_set(bf, 0, 50);
    nlme_bitfield_set(bf, 100, 150);
    nlme_bitfield_unset(bf, 50, 100);
    nlme_bitfield_unset(bf, 150, 200);
    nlme_bitfield_set(bf, 200, 250);
    nlme_bitfield_set(bf, 300, 350);
    nlme_bitfield_unset(bf, 250, 300);
    nlme_bitfield_unset(bf, 350, 400);
    nlme_bitfield_set(bf, 155, 156);
    nlme_bitfield_unset(bf, 305, 306);
    nlme_bitfield_print(bf, 0, 400);

    range.end = 0;
    while (nlme_bitfield_next_range(bf, 400, &range)) {
        printf("[%d, %d)\n", range.beg, range.end);
    }
    for (i = 0; i < 10000000; i++) {
        range.end = 0;
        while (nlme_bitfield_next_range(bf, 400, &range)) {
        }
    }
    nlme_bitfield_free(bf);
    bf = nlme_bitfield_alloc(len);
    nlme_bitfield_unset(bf, 0, len-1);
    nlme_bitfield_set(bf, len-1, len);
    range.end = 0;
    nlme_bitfield_next_range(bf, len, &range);
    printf("[%d, %d)\n", range.beg, range.end);
    return 0;
}
#endif

#if (CHAR_BIT != 8)
#error "System does not have exactly 8 bits in char types"
#endif

/** 
 * Allocate a new nlme_bitfield object.
 * 
 * @param len number of bits to be stored
 * 
 * @return newly allocated nlme_bitfield object
 */
nlme_bitfield *
nlme_bitfield_alloc(nlme_word len)
{
    nlme_bitfield *ans = (nlme_bitfield *)
        Calloc(sizeof(nlme_bitfield)+
               sizeof(int)*((len+WORD_BITS-1)/WORD_BITS), char);
    if (ans == NULL)
        error("could not allocate bitfield");
    ans->len = len;
    ans->words = (nlme_word *) (ans + 1);
    return ans;
}

/** 
 * Free an allocated nlme_bitfield object
 * 
 * @param x object to be freed
 */
void
nlme_bitfield_free(nlme_bitfield *x)
{
    Free(x);
}

#define ALL_ONE (~((nlme_word) 0))

/** 
 * Toggle a range of indices in a bitfield
 * 
 * @param x nlme_bitfield object
 * @param start start of the range
 * @param end start of the range
 */
void
nlme_bitfield_toggle(nlme_bitfield *x, nlme_word beg,
                     nlme_word end)
{
    nlme_word i;
    nlme_word beg1, end1;

    if (x == NULL || beg == end)
        return;

    if (x->len <= beg || x->len < end || end < beg)
        error("Invalid index to bitfield. beg: %d, end: %d, length: %d", beg, end, x->len);

    beg1 = beg/WORD_BITS;
    end1 = (end-1)/WORD_BITS;
    beg %= WORD_BITS;
    end %= WORD_BITS;
    if (beg1 == end1) {
        x->words[beg1++] ^=
            (ALL_ONE << (WORD_BITS+beg-end)) >> (WORD_BITS-end);
    } else {
        if (beg > 0)
            x->words[beg1++] ^= ALL_ONE << (beg);

        if (end > 0)
            x->words[end1] ^= ALL_ONE >> (WORD_BITS-end);
        else end1++;

        for (i = beg1; i < end1; i++)
            x->words[i] ^= ALL_ONE;
    }
}

/** 
 * Set a range of indices in a bitfield
 * 
 * @param x nlme_bitfield object
 * @param start start of the range
 * @param end start of the range
 */
void
nlme_bitfield_set(nlme_bitfield *x, nlme_word beg,
                  nlme_word end)
{
    nlme_word beg1;
    nlme_word end1;

    if (x == NULL || beg == end)
        return;

    if (x->len <= beg || x->len < end || end < beg)
        error("Invalid index to bitfield. beg: %d, end: %d, length: %d", beg, end, x->len);

    beg1 = beg/WORD_BITS;
    end1 = (end-1)/WORD_BITS;
    beg %= WORD_BITS;
    end %= WORD_BITS;
    if (beg1 == end1) {
        x->words[beg1++] |=
            (ALL_ONE << (WORD_BITS+beg-end)) >> (WORD_BITS-end);
    } else {
        if (beg > 0)
            x->words[beg1++] |= ALL_ONE << (beg);

        if (end > 0)
            x->words[end1] |= ALL_ONE >> (WORD_BITS-end);
        else end1++;

        memset(&(x->words[beg1]), 0xff, (end1-beg1)*sizeof(nlme_word));
    }
}

/** 
 * Unset a range of indices in a bitfield
 * 
 * @param x nlme_bitfield object
 * @param start start of the range
 * @param end start of the range
 */
void
nlme_bitfield_unset(nlme_bitfield *x, nlme_word beg,
                    nlme_word end)
{
    nlme_word beg1;
    nlme_word end1;

    if (x == NULL || beg == end)
        return;

    if (x->len <= beg || x->len < end || end < beg)
        error("Invalid index to bitfield. beg: %d, end: %d, length: %d", beg, end, x->len);

    beg1 = beg/WORD_BITS;
    end1 = (end-1)/WORD_BITS;
    beg %= WORD_BITS;
    end %= WORD_BITS;
    if (beg1 == end1) {
        x->words[beg1++] &=
            ~(ALL_ONE << (WORD_BITS+beg-end) >> (WORD_BITS-end));
    } else {
        if (beg > 0)
            x->words[beg1++] &= ALL_ONE >> (WORD_BITS-beg);

        if (end > 0)
            x->words[end1] &= ALL_ONE << end;
        else end1++;

        memset(&(x->words[beg1]), 0, (end1-beg1)*sizeof(nlme_word));
    }
}

/** 
 * Return the base-2 integer logarithm of the word.
 *
 * We assume that the word has exactly one non-zero bit.
 * 
 * @param word an unsigned integer value
 * 
 * @return base-2 logarithm of the word.
 */
static nlme_word
nlme_word_log2(nlme_word x)
{
    nlme_word m;
#if defined(__GNUC__) && defined(i386)
    __asm__("bsrl %1,%0\n\t"
            : "=r" (m)
            : "g"  (x));
    return m;
#else
#if 1
    nlme_word n;
    if (sizeof(nlme_word) > 4) { /* 64 bit word */
        n = ((x >> (sizeof(nlme_word) > 4)?32:0) != 0) << 5;
        x >>= n;
    } else n = 0;
    m = ((x >> 16) != 0) << 4;
    n += m;
    x >>= m;
    m = ((x >> 8) != 0) << 3;
    x >>= m;
    n += m;
    m = ((x >> 4) != 0) << 2;
    x >>= m;
    return (n + m + ((804>>x)&(x-1)));
#else
    signed char first_set_bit[] = {
       -1, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
        4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0
    };

    int i = 0;
    unsigned int byte = x & 0xff;

    if (byte)
        return first_set_bit[byte];
    byte = x & (0xff << 8);
    if (byte)
        return first_set_bit[byte >> 8]+CHAR_BIT;
    byte = x & (0xff << 16);
    if (byte)
        return first_set_bit[byte >> 16]+2*CHAR_BIT;
    byte = x & (0xff << 24);
    if (byte)
        return first_set_bit[byte >> 24]+3*CHAR_BIT;
    if (sizeof(int) > 4) { /* should be optimized away for machines
                            * for which this is false */
        x >>= (sizeof(int)>4)?32:0;
        return 4*CHAR_BIT+nlme_fsbit(word);
    }
#endif
#endif
}

/** 
 * 
 * 
 * @param x nlme_bitfield object
 * @param end range should be below this
 * @param range on input the last range returned, on output the next range
 * 
 * @return 1 if a new range was found, 0 otherwise.
 */
int
nlme_bitfield_next_range(const nlme_bitfield *x, nlme_word end,
                         nlme_range *range)
{
    nlme_word beg;

    if (range == NULL)
        error("range pointer must be non NULL");

    beg = range->end;
    if (x == NULL) {
        if (beg >= end)
            goto no_range_found;
        range->beg = beg;
        range->end = end;
        return 1;
    }

    if (end > x->len)
        end = x->len;
    if (beg < end) {
        nlme_word i = beg/WORD_BITS;
        nlme_word end1 = 1+(end-1)/WORD_BITS;
        nlme_word word = x->words[i++]>>(beg%WORD_BITS);
        nlme_word tmp;

        if (!word) {
            for (; i < end1; i++) {
                word = x->words[i];
                if (word) {
                    beg = i*WORD_BITS;
                    i++;
                    goto have_begining;
                }
            }
            goto no_range_found;
        }
    have_begining:
        tmp = word;
        word &= -word;
        range->beg = beg + nlme_word_log2(word);
        if (end <= range->beg)
            goto no_range_found;
        word += tmp;
        if (!word) {
            for (; i < end1; i++) {
                word = x->words[i] + 1;
                if (word) {
                    beg = i*WORD_BITS;
                    goto have_ending;
                }
            }
            range->end = end;
            return 1;
        }
    have_ending:
        word &= -word;
        word = beg + nlme_word_log2(word);
        range->end = min(end, word);
        return 1;
    }
no_range_found:
    range->beg = range->end = end;
    return 0;
}
