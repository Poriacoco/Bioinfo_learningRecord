#include<stdio.h>
#include<unistd.h>

int main()
{
    int cnt=0;
    while(1)
    {
        int ret=fork();
        if(ret<0){
            printf("error,%d\n",cnt);
            break;
        }
        if(ret==0)
            while(1) sleep(1);
    cnt++;
    }
     return 0;
}
