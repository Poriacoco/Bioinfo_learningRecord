#include "process.h"
char style[StyNum] = {'-','#','+'};

void ProcessOn()//函数定义

{
    int cnt =0;
    char bar[NUM];
    memset(bar,'\0',sizeof(bar));

    while(cnt<=100)
    {
       printf("[%-100s][%3d%%]\r",bar,cnt);
       //缓冲区刷新数据
       fflush(stdout);
       bar[cnt++]=style[N];//N会在makefile文件中定义
       usleep(50000);
    }
    printf("\n");
    return ;

}










