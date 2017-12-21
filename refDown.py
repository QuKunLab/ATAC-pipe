import os
os.popen("findMotifs.pl 1> findMotifs.pl.log 2> findMotifs.pl.err")
errfile=open("findMotifs.pl.err",'r')
lines=errfile.readlines()
cmd=''
for line in lines:
    if line.split('"')[-1]==' to install promoter set NNN\n':
        cmd=cmd+(line.split('"')[1])[:-3]
        print cmd
cmd1=cmd+"human-p"
cmd2=cmd+"human-o"
cmd3=cmd+"mouse-p"
cmd4=cmd+"mouse-o"
cmd5=cmd+"hg19"
cmd6=cmd+"mm9"
os.popen(cmd1)
os.popen(cmd2)
os.popen(cmd3)
os.popen(cmd4)
os.popen(cmd5)
os.popen(cmd6)
os.popen("rm findMotifs.pl.log")
os.popen("rm findMotifs.pl.err")

os.popen("sh ./Code/chr.sh")

