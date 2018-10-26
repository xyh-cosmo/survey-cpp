#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>

#define MAXLINE 2046

int main(int argc, char *argv[]) {
    FILE* file = fopen("", "r");
    char line[MAXLINE];
    while (fgets(line, MAXLINE, file) != NULL ) {

    }
}
