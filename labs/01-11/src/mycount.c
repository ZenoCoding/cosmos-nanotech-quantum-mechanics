#include <stdio.h>
#include <math.h>

int main(void){
    int x;
    printf("How high? ");
    scanf("%d", &x);
    for(int i = 0; i <= x; i++){
        printf("%d ", i);
    }
    return 0;
}