#include<stdio.h>
#include<unistd.h>

int main()
{
    pid_t id =fork();
    if(id < 0)
    {
        perror("fork");
        return 1;
    }
    else if(id == 0)
    {
        //child
        while(1)
        {
            printf("child,pid: %d,ppid: %d\n",getpid(),getppid());
            sleep(1);
        }
    }
    else  
    {
        //parent
        while(2)
         {
            printf("parent,pid: %d,ppid: %d\n",getpid(),getppid());
            sleep(1);
        }
    }

}

