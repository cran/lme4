#ifndef BITFIELD_H
#define BITFIELD_H

typedef unsigned int nlme_word;

typedef struct {
    nlme_word beg;
    nlme_word end;
} nlme_range;

struct nlme_bitfield; /* opaque type for bitfields */
typedef struct nlme_bitfield nlme_bitfield;

nlme_bitfield *nlme_bitfield_alloc(nlme_word len);
void nlme_bitfield_free(nlme_bitfield *x);
void nlme_bitfield_print(nlme_bitfield *x, nlme_word start,
                         nlme_word end);
void nlme_bitfield_toggle(nlme_bitfield *x, nlme_word start,
                          nlme_word end);
void nlme_bitfield_set(nlme_bitfield *x, nlme_word start,
                       nlme_word end);
void nlme_bitfield_unset(nlme_bitfield *x, nlme_word start,
                         nlme_word end);
/**
 * This can be called like this:
 *    nlme_range range;
 *    int end = SOME_VALUE;
 *    while (nlme_bitfield_next_range(x, end, &range)) {
 *        // do something with the range
 *    }
 */

int nlme_bitfield_next_range(const nlme_bitfield *x,
                             nlme_word end, nlme_range* range);
#endif
