#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main()
{
 pid_t id = fork();
if(id > 0)
 { //parent
    printf("parent[%d] is sleeping...\n", getpid());
    sleep(5);
 }
 else
 {
     //child
    printf("child[%d] is begin Z...\n", getpid());
    sleep(1);
    exit(EXIT_SUCCESS);
 }
 return 0;
 }
