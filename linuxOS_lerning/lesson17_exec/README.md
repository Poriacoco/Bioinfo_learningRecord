### myexec
- myexec* 包含5种常见的程序调用接口
- myexec.c中使用了最基本的execl 实现了利用fork()创建子进程,实现ls -a -l的调用

### myshell
- myshell 目录中利用execvp程序调用接口实现了一个简单的shell
- 完成了ls catch touch top 等基础命令的实现
