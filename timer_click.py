# -*- coding:utf-8 -*-
import time
from pymouse import  PyMouse

if __name__=='__main__':
    index=0
    sleeps=raw_input("click interval(s):")
    while(not sleeps.isdigit()):
        sleeps=raw_input("click interval(s):")
    sleeps=int(sleeps)    
    click_num=raw_input("click_num:")
    while(not click_num.isdigit()):
        click_num=raw_input("click_num:")
    click_num=int(click_num)
    points=[]
    key=""
    while(index<click_num and key!='q' and key!='Q'):
        key=raw_input('press w get click %s pos,or q/Q for quit:'%(index+1))
        while(key=='w' or key=='W'):
            m=PyMouse()
            x,y=m.position()
            points.append((x,y))
            index+=1
            key=""
    click_num=index
    index=0
    while(True):
        m=PyMouse()
        x,y=points[index]
        print time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())),'click at(',x,',',y,')'
        m.move(x,y)
        m.click(x,y)
        index+=1
        if(index>=click_num):
            break
        time.sleep(sleeps)

    


 
        