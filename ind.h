//
// Created by fede on 21/01/19.
//

#ifndef LEVENSHTEIN_IND_H
#define LEVENSHTEIN_IND_H

typedef struct ind{
    ind(){
        i = 0; j = 0;
    }
    ind(int i, int j): i(i), j(j){};
    /*ind(ind& ind1) {
        i = ind1.i;
        j = ind1.j;
    };*/
    int i, j;
} ind;

#endif //LEVENSHTEIN_IND_H
