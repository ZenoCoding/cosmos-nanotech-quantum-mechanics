#include <stdio.h>
#include <math.h>

#define ll long long

int main(void){
    ll ct;
    printf("How many do I sum? ");
    scanf("%llu", &ct);
    ll total = 0;
    for(ll i = 0; i < ct; i++){
        total += i + 1;
    }
    printf("Answer: %llu\n", total);
}