#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <assert.h>

#define NUM 1024
#define OPT_NUM 64

char lineCommand[NUM];
char *myargv[OPT_NUM]; //指针数组

int main()
{

while(1)//死循环,保证命令行可以不断进行输入
    {
    //输入提示符
    printf("user@host PATH# ");
    fflush(stdout);//刷新缓冲区

    //获取用户输入
    char *s = fgets(lineCommand,sizeof(lineCommand)-1,stdin);
   assert(s!=NULL);
   (void)s;
   //输入命令的时候最后一个输入一定是回车键入,故需要处理最后一个\n
   lineCommand[strlen(lineCommand)-1]=0;
   /*printf("test:%s\n",lineCommand);*/ //测试代码

   //字符串切割 "ls -a -l "->"ls" "-a" "-l"
   myargv[0]=strtok(lineCommand," ");//按照空格切割第一个,ls 
    //如果后续没有子字符串,strtok最后会返回NULL,故myargv[end]=NULL;
    int i=1;
    //对字符串反复切割并将其存入数组myargv中,且保证myargv[end]=NULL;
    while(myargv[i++]=strtok(NULL," "));
    //对字符串做继续切割,strtok传参值为NULL
    
    //测试是否成功  (条件编译) *Makefile 添加-DEBUG
    #ifdef DEBUG 
    for(int i=0;myargv[i];i++)
        printf("myargv[%d]: %s\n",i,myargv[i]);
    #endif

    //fork创建子进程,execvp执行函数调用
    pid_t id=fork();
    assert(id!= -1);
    if(id ==0)
    {
        execvp(myargv[0],myargv);//'v'表示vector,命令是myargv数组里的一串参数;'p'表示PATH,函数会自己寻找PATH环境变量,无需键入
        //myargv[end]==NULL,自动传入execvp中,运行成功则进程结束,运行失败会有返回值,需要exit(1)执行退出和退出码
        exit(1);//执行失败必定返回
    }
    waitpid(id,NULL,0);
    }
}








