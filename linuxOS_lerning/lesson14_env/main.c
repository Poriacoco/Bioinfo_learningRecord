#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define USER "USER"
#define MY_NEV "myval"
#define MYPWD "PWD"

// ls -a -l 选项个数-int argc//
// 选项输入"-a""-l"的地址在*argv[]的指针数组里
int main(int argc,char *argv[])
{
    if(argc!=2)//选项个数不满足条件
    {
        printf("Usage:\n\t%s[-a/-l/-al]\n",argv[0]);
        //argv[0]是指令,argv[1...2....3...]是命令行参数
        return 0;
    }

    if(strcmp("-a",argv[1])==0)
    {
        printf("function A");
    }
    if(strcmp("-l",argv[1])==0)
    {
        printf("function L");
    }
    if(strcmp("-al",argv[1])==0)
    {
        printf("function AL");
    }
        return 0;
}
