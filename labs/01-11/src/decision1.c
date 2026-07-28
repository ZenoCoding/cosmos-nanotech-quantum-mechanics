/*number comparer*/
#include <stdio.h>
#include <math.h>

int main(){
    float x, y;
    printf("Enter x and y: ");
    scanf("%f %f", &x, &y);
    if(x > y)
    printf("first number bigger");
    else printf("second bigger");
    return 0;
}