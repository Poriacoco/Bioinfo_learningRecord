########################################################
#                  六十三 安装操作系统                     #
########################################################
#1 CentOS
https://www.centos.org/

#2 Ubuntu
http://www.ubuntu.com/ 

#3 腾讯云
#链接地址：https://curl.qcloud.com/gm6m0QoY


########################################################
#                   六十四 磁盘配置                     #
########################################################
#=========================
#       1 分区格式化      #
#=========================
#创建一个文件夹，用户挂载磁盘
mkdir /ifs1
#进行格式化
fdisk -l
parted /dev/vdb  #交互界面 mklabel gpt ，quit
mkfs.xfs -f /dev/vdb 
#挂载磁盘
mount /dev/vdb /ifs1
#设置自动挂载，将下面信息追加写入/etc/fstab文件中
echo "/dev/vdb /ifs1                       xfs     defaults,uquota        0 0" >>/etc/fstab


#=========================
#       2 磁盘管理        #
#=========================
mkdir -p /ifs1/Software/src  /ifs1/Software/biosoft  /ifs1/Software/bin
mkdir User
mkdir Database
