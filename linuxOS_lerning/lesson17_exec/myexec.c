#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/wait.h>
int test()
{
    //c语言程序替换
    //流程 .c~exe~load~process~running
    
    printf("process is running...\n");
//load
    execl("/usr/bin/ls","ls","-a","-l",NULL);
//int execl(const char*path,const char *arg,...); 所有的exec函数均以NULL结尾
 
    printf("process running done...\n");

    return 0;
        
}

int main()
{
    //子进程执行程序替换
    
    printf("process is running...\n");
    pid_t id=fork();
    assert(id!=-1);

    if(id==0)
    {
        sleep(1);
        execl("/usr/bin/ls","ls","-a","-l",NULL);
        exit(1);/*execl一旦退出必定是异常退出,执行成功不退出*/
    }
    int status=0;
    pid_t ret=waitpid(id,&status,0);
    if(ret>0)
    printf("wait success:exit code:%d,sig:%d\n",(status>>8)&0xFF,status & 0x7F);
}

