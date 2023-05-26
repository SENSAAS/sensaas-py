#!/usr/bin/python3.7

#Author: Dominique DOUGUET, Lucas GRANDMOUGIN, Mohamed-Akram MASROUR and Hamza BIYUZAN

import sys
import os
import math
import re
import numpy as np
#from scipy.spatial import KDTree
#from biopandas.pdb import PandasPdb

##########################################
def molsurface(filetype,molfile,TemplateDots,nbDots,printdots):

    #label1 {H, Cl, Br, I} white/grey 0.9 0.9 0.9
    #label2 {O, N, S, F} red 1 0 0
    #label3 {C, P, B} green 0 1 0
    #label4 {others} blue 0 0 1

    tabR= {'C':'%.2f' % 1.70, 'O':1.52, 'N':1.55, 'S':1.80, 'P':1.80, 'B':1.72, 'Br':1.85, 'Cl':1.75, 'I':1.98, 'F':1.47, 'H':'%.2f' % 1.20, 'Hp':'%.2f' % 1.10, 'X':'%.2f' % 1.10}
    #no label for 'X' Dummy atoms : no associated dots are saved = hole in point cloud
    label= {'C':3, 'P':3, 'B':3, 'O':2, 'N':2, 'S':2, 'F':2, 'Hp':2, 'H':1, 'Cl':1, 'Br':1, 'I':1, 'X':1}
    rgb= np.array([[0, 0, 0], [0.9, 0.9, 0.9], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    espace5=' '
    espace6='       '
    fichier2D=0
    verbose=0
    if(printdots==1):
        verbose=1
    #either uses KDTree (1) or distance calculation when atoms radius differ (0) (statistically it takes longer)
    usekdtree=0
    #if = 1 it reads pdb using biopandas (statistically it takes longer)
    usepandaspdb=0
    
    #open mol file 
    if(filetype=="sdf" or (filetype=="pdb" and usepandaspdb==0)):
        filemol=open(molfile,'r')
        getstr=filemol.read().split('\n')
        filemol.close()

################# open and read sdf

    if(filetype=="sdf"):
        tabLignesSdf=[]
        compt=3
        while(compt < len(getstr)):
            tabLignesSdf.append(re.split('\s+', getstr[compt].strip()))
            compt=compt+1
        testspace=[]
        testspace.append(re.split('', getstr[3]))
        if(len(tabLignesSdf[0][0]) > 2):#attached
            #" 11222" or "111 22" or "111222"
            if(testspace[0][1]==' '):#" 11222"
                #keep the following order:
                tabLignesSdf[0][1]=tabLignesSdf[0][0][2:]
                tabLignesSdf[0][0]=tabLignesSdf[0][0][0:2]
            elif(testspace[0][4]!=' '):
                #keep the following order:
                tabLignesSdf[0][1]=tabLignesSdf[0][0][3:]
                tabLignesSdf[0][0]=tabLignesSdf[0][0][:3]
            nbatomes=int(tabLignesSdf[0][0])
            nbLiaisons=int(tabLignesSdf[0][1])
        else:#" 11 22" or "11 222"
            nbatomes=int(tabLignesSdf[0][0])
            nbLiaisons=int(tabLignesSdf[0][1])

        #Extract coordinates, atom type and set radius
        compt=1
        getx=[]
        getx.append('')
        gety=[]
        gety.append('')
        getz=[]
        getz.append('')
        getA=[]
        getA.append('')
        getRayon=[]
        getRayon.append('')
        getColor=[]
        getColor.append('')
        arr_xyz=np.empty(shape=[nbatomes,4], dtype='float64')
        while (compt <= nbatomes): 
            arr_xyz[compt-1,0]=float(tabLignesSdf[compt][0])
            arr_xyz[compt-1,1]=float(tabLignesSdf[compt][1])
            arr_xyz[compt-1,2]=float(tabLignesSdf[compt][2])
            getx.append(arr_xyz[compt-1,0])
            gety.append(arr_xyz[compt-1,1])
            getz.append(arr_xyz[compt-1,2])
            if (arr_xyz[compt-1,2] == 0):
                fichier2D=fichier2D+1
            getA.append(tabLignesSdf[compt][3])
            if(getA[compt] not in tabR):
                print("Warning: atom %s set as C because it is not the tab (unusual in medchem)" % getA[compt])
                getA[compt]='C'
            arr_xyz[compt-1,3]=float(tabR[getA[compt]])
            getRayon.append(float(tabR[getA[compt]]))
            getColor.append(label[getA[compt]])
            compt=compt+1
        
        #Search polar H
        maxl=nbatomes+nbLiaisons+4
        compt=nbatomes+4
        while (compt < maxl):
            atom1=int(getstr[compt][:3])
            atom2=int(getstr[compt][3:6])
            if (atom2 > nbatomes or atom2 < 0 or atom1 > nbatomes or atom1  < 0):
                print("invalid atom number %6d or %6d" % (atom1,atom2))
                quit()
            if (((getA[atom1] == 'O') and (getA[atom2] == 'H')) or (getA[atom1] == 'H') and (getA[atom2] == 'O')):
                if(getA[atom1]=='H'):
                    getRayon[atom1]=float(tabR['Hp'])
                    getColor[atom1]=label['Hp']
                    arr_xyz[atom1-1,3]=float(tabR['Hp'])
                elif(getA[atom2]=='H'):
                    getRayon[atom2]=float(tabR['Hp'])
                    getColor[atom2]=label['Hp']
                    arr_xyz[atom2-1,3]=float(tabR['Hp'])
            if (((getA[atom1] == 'N') and (getA[atom2] == 'H')) or (getA[atom1] == 'H') and (getA[atom2] == 'N')):
                if(getA[atom1]=='H'):
                    getRayon[atom1]=float(tabR['Hp'])
                    getColor[atom1]=label['Hp']
                    arr_xyz[atom1-1,3]=float(tabR['Hp'])
                elif(getA[atom2]=='H'):
                    getRayon[atom2]=float(tabR['Hp'])
                    getColor[atom2]=label['Hp']
                    arr_xyz[atom2-1,3]=float(tabR['Hp'])
            compt=compt+1 

############# open and read pdb

    elif(filetype=="pdb" and usepandaspdb==0):
        tabLignesPdb=[]
        tabLignesPdb.append('')
        compt=1
        while (compt < len(getstr)):
            tabLignesPdb.append(re.split('\s+', getstr[compt].strip()))
            compt=compt+1
        compt=1
        nbatomes=0
        getx=[]
        getx.append('')
        gety=[]
        gety.append('')
        getz=[]
        getz.append('')
        getA=[]
        getA.append('')
        getRayon=[]
        getRayon.append('')
        getColor=[]
        getColor.append('')
        while (compt < len(tabLignesPdb)):
            if (tabLignesPdb[compt][0] == 'HETATM' or tabLignesPdb[compt][0] == 'ATOM'):
                xAtome=float(tabLignesPdb[compt][5])
                yAtome=float(tabLignesPdb[compt][6])
                zAtome=float(tabLignesPdb[compt][7])
                getx.append(xAtome)
                gety.append(yAtome)
                getz.append(zAtome)
                if (float(zAtome) == 0):
                    fichier2D=fichier2D+1
                getA.append(tabLignesPdb[compt][2])
                if(getA[compt] not in tabR):
                    print("Warning: atom %s set as C because it is not the tab (unusual in medchem)" % getA[compt])
                    getA[compt]='C'
                getRayon.append(float(tabR[getA[compt]]))
                getColor.append(label[getA[compt]])
                nbatomes=nbatomes+1
            compt=compt+1

        compt=1
        arr_xyz=np.empty(shape=[nbatomes,4], dtype='float64')
        while (compt <= nbatomes):
             arr_xyz[compt-1,0]=getx[compt]
             arr_xyz[compt-1,1]=gety[compt]
             arr_xyz[compt-1,2]=getz[compt]
             arr_xyz[compt-1,3]=getRayon[compt]
             compt=compt+1

############# or uses BioPandas

    elif(filetype=="pdb" and usepandaspdb==1):
        ppdb = PandasPdb()
        ppdb.read_pdb(molfile)
        testhetatm =  ppdb.df['HETATM'].empty
        testatom = ppdb.df['ATOM'].empty
        if(not testhetatm and testatom):
        #HETATM only
            datah = ppdb.df['HETATM']
            #xyzh = datah[['x_coord','y_coord','z_coord','element_symbol']]
            xyzh = datah[['x_coord','y_coord','z_coord','atom_name']]
            arr_xyzb = np.array(xyzh)
        elif(not testatom and testhetatm==1):
        #ATOM only
            data1 = ppdb.df['ATOM']
            #warning: element_symbol is not always indicated !
            xyz1 = data1[['x_coord','y_coord','z_coord','element_symbol']]
            arr_xyzb = np.array(xyz1)
        elif(not testhetatm and not testatom):
        #HETATM and ATOM
            datah = ppdb.df['HETATM']
            xyzh = datah[['x_coord','y_coord','z_coord','atom_name']]
            arr_xyzh = np.array(xyzh)
            data1 = ppdb.df['ATOM']
            xyz1 = data1[['x_coord','y_coord','z_coord','element_symbol']]
            arr_xyz1 = np.array(xyz1)
        #here choose to keep hetatm only or both:
            #arr_xyzb = np.concatenate((arr_xyz1,arr_xyzh), axis=0)
            arr_xyzb = np.copy(arr_xyzh)
        
        widthxyz = len(arr_xyzb)
        compt=1
        nbatomes=0
        getx=[]
        getx.append('')
        gety=[]
        gety.append('')
        getz=[]
        getz.append('')
        getA=[]
        getA.append('')
        getRayon=[]
        getRayon.append('')
        getColor=[]
        getColor.append('')
        while (nbatomes < widthxyz):
            xAtome=float(arr_xyzb[nbatomes,0])
            yAtome=float(arr_xyzb[nbatomes,1])
            zAtome=float(arr_xyzb[nbatomes,2])
            getx.append(xAtome)
            gety.append(yAtome)
            getz.append(zAtome)
            if (float(zAtome) == 0):
                fichier2D=fichier2D+1
            getA.append(str(arr_xyzb[nbatomes,3]))
            if(getA[compt] not in tabR):
                print("Warning: atom %s set as C because it is not the tab (unusual in medchem)" % getA[compt])
                getA[compt]='C'
            getRayon.append(float(tabR[getA[compt]]))
            getColor.append(label[getA[compt]])
            nbatomes=nbatomes+1
            compt=compt+1

        compt=1
        arr_xyz=np.empty(shape=[nbatomes,4], dtype='float64')
        while (compt <= nbatomes):
             arr_xyz[compt-1,0]=getx[compt]
             arr_xyz[compt-1,1]=gety[compt]
             arr_xyz[compt-1,2]=getz[compt]
             arr_xyz[compt-1,3]=getRayon[compt]
             compt=compt+1

    #Search polar in pdb
    if(filetype=="pdb"):
        compt=1
        while (compt <= nbatomes):
            if (getA[compt] == 'H'):
                compt2=1
                while(compt2 <= nbatomes):
                    if (getA[compt2] == 'N' or getA[compt2] == 'O'):
                        distHp= math.sqrt((getx[compt] - getx[compt2])**2 + (gety[compt] - gety[compt2])**2 + (getz[compt] - getz[compt2])**2)
                        if (distHp <= 1.2):
                           getRayon[compt]=float(tabR['Hp'])
                           getColor[compt]=label['Hp']
                    compt2=compt2+1
            compt=compt+1

########## end pdb

########## checks

    if (fichier2D==int(nbatomes)):
        print("Warning: mol file in 2D; SenSaaS needs 3D coordinates to work properly")

########## box size

    tmin=np.min(arr_xyz, axis=0)
    tmax=np.max(arr_xyz, axis=0)
    minx=tmin[0]
    maxx=tmax[0]
    miny=tmin[1]
    maxy=tmax[1]
    minz=tmin[2]
    maxz=tmax[2]
    rmax=tmax[3]
        
########### end read molfile

    #initialize KDtree for speeding neighbour atom searches
    if(usekdtree == 1):
        atoms = np.zeros((nbatomes,3))
        atoms[:, 0] = np.array(getx[1:])
        atoms[:, 1] = np.array(gety[1:])
        atoms[:, 2] = np.array(getz[1:])
        tree = KDTree(atoms)

    usetemplate=1
    if(usetemplate):
        #check if surface dots are in a cell = atomxyz +/- radius in each cell of the main grid
        deuxrmax=2*rmax
        boxx=maxx-minx+deuxrmax
        boxy=maxy-miny+deuxrmax
        boxz=maxz-minz+deuxrmax
        if (boxx >= boxy and boxx >= boxz):
            sizebox=boxx
        elif (boxy >= boxx and boxy >= boxz):
            sizebox=boxy
        else:
            sizebox=boxz
        nbbox=int(sizebox/deuxrmax)+1
        nbcells=nbbox*nbbox*nbbox
        lengthbox=deuxrmax
        bl=1
        indicebl=0
        maxblist=nbcells*nbatomes
        #boxi, atomi
        blistbox=np.empty(shape=[maxblist,2], dtype='int')
        xlow=minx-rmax
        xhigh=xlow+lengthbox
        for bi in range(1,nbbox+1):
            ylow=miny-rmax
            yhigh=ylow+lengthbox
            for bj in range(1,nbbox+1):
                zlow=minz-rmax
                zhigh=zlow+lengthbox
                for bk in range(1,nbbox+1):
                    compt=1
                    while (compt <= nbatomes):
                        if ( ((getx[compt]+getRayon[compt]) <= xhigh and (getx[compt]+getRayon[compt]) >= xlow) or ((getx[compt]-getRayon[compt]) <= xhigh and (getx[compt]-getRayon[compt]) >= xlow) ):
                            if ( ((gety[compt]+getRayon[compt]) <= yhigh and (gety[compt]+getRayon[compt]) >= ylow) or ((gety[compt]-getRayon[compt]) <= yhigh and (gety[compt]-getRayon[compt]) >= ylow) ):
                                if ( ((getz[compt]+getRayon[compt]) <= zhigh and (getz[compt]+getRayon[compt]) >= zlow) or ((getz[compt]-getRayon[compt]) <= zhigh and (getz[compt]-getRayon[compt]) >= zlow) ):
                                    blistbox[indicebl,:]=[bl, compt]
                                    indicebl=indicebl+1
                        compt=compt+1
                    bl=bl+1
                    zlow=zhigh
                    zhigh=zhigh+lengthbox
                ylow=yhigh
                yhigh=yhigh+lengthbox
            xlow=xhigh
            xhigh=xhigh+lengthbox

        versionbool=0
        if(versionbool==1):
            #Matrix: set True if sphere (part of) in same box
            blistboxb=np.empty(shape=[nbatomes+1,nbatomes+1], dtype='bool')
            for bl in range(0,indicebl):
                for bl2 in range(0,indicebl):
                    if(bl!=bl2 and blistbox[bl,0]==blistbox[bl2,0]):
                        blistboxb[blistbox[bl,1],blistbox[bl2,1]]=True
                        blistboxb[blistbox[bl2,1],blistbox[bl,1]]=True

        #maximum:
        nbtotalDots=nbDots*nbatomes
        getDots=np.empty(shape=[nbtotalDots,3], dtype='float64')
        getrgb=np.empty(shape=[nbtotalDots,3], dtype='float64')
        comptDots=0
        labeltot=[]
        label1=[]
        label2=[]
        label3=[]
        label4=[]
        #remove masked dots
        compt=1
        while (compt <= nbatomes):
            #register products for atom compt
            prod=np.zeros(shape=[(nbatomes+1),5], dtype='float64')
            prodi=0
            if(versionbool==1):
                for compt2 in range(1,nbatomes):
                    if(blistboxb[compt,compt2]):
                        ajix=getx[compt2]-getx[compt]
                        ajixsq=ajix*ajix
                        ajiy=gety[compt2]-gety[compt]
                        ajiysq=ajiy*ajiy
                        ajiz=getz[compt2]-getz[compt]
                        ajizsq=ajiz*ajiz
                        dij=ajixsq+ajiysq+ajizsq
                        ri2=getRayon[compt]*getRayon[compt]
                        rj2=getRayon[compt2]*getRayon[compt2]
                        limit=(dij+ri2-rj2)/(2*getRayon[compt])
                        prod[prodi,:]=[ajix, ajiy, ajiz, limit, compt2]
                        prodi=prodi+1
            else:
                vu=np.zeros(shape=[(nbatomes+1)], dtype='int32')
                for bl in range(0,indicebl):
                    if(blistbox[bl,1]==compt):
                        for bl2 in range(0,indicebl):
                            compt2=int(blistbox[bl2,1])
                            if ( blistbox[bl,0]==blistbox[bl2,0] and bl!=bl2 and vu[compt2]==0 ):
                                ajix=getx[compt2]-getx[compt]
                                ajixsq=ajix*ajix
                                ajiy=gety[compt2]-gety[compt]
                                ajiysq=ajiy*ajiy
                                ajiz=getz[compt2]-getz[compt]
                                ajizsq=ajiz*ajiz
                                dij=ajixsq+ajiysq+ajizsq
                                ri2=getRayon[compt]*getRayon[compt]
                                rj2=getRayon[compt2]*getRayon[compt2]
                                limit=(dij+ri2-rj2)/(2*getRayon[compt])
                                prod[prodi,:]=[ajix, ajiy, ajiz, limit, compt2]
                                prodi=prodi+1
                                vu[compt2]=1

            #dotx, doty, dotz, kept
            tabxyz=np.empty(shape=[nbDots], dtype='float64')
            for tabi in range(0,nbDots):
                tabxyz[tabi]=-1
            for tabi in range(0,nbDots):
                if (tabxyz[tabi]==-1):
                    flag=0
                    for j in range(0,prodi):
                        projdotx=TemplateDots[tabi,0]*prod[j,0]
                        projdoty=TemplateDots[tabi,1]*prod[j,1]
                        projdotz=TemplateDots[tabi,2]*prod[j,2]
                        projdot=projdotx+projdoty+projdotz
                        if (projdot > prod[j,3]):
                            tabxyz[tabi]=0
                            flag=1
                            for k in range(tabi+1,nbDots):
                                projdotx=TemplateDots[k,0]*prod[j,0]
                                projdoty=TemplateDots[k,1]*prod[j,1]
                                projdotz=TemplateDots[k,2]*prod[j,2]
                                projdot=projdotx+projdoty+projdotz
                                if (projdot > prod[j,3]):
                                    tabxyz[k]=0
                                else:
                                    break    
                        if (flag==1):
                            break
                    if (flag==0):
                        getDots[comptDots,0]='%.2f' % (getRayon[compt] * TemplateDots[tabi,0] + getx[compt])
                        getDots[comptDots,1]='%.2f' % (getRayon[compt] * TemplateDots[tabi,1] + gety[compt])
                        getDots[comptDots,2]='%.2f' % (getRayon[compt] * TemplateDots[tabi,2] + getz[compt])
                        #set color
                        mi=compt
                        #Search KDtree for the closest atom
                        if(usekdtree == 1):
                            qi = tree.query(getDots[comptDots,:])
                            mi = qi[1] + 1
                        else:
                            closestR=getRayon[compt]
                            compt2=1
                            while (compt2 <= nbatomes):
                                #case radius atom compt2 (eg H) < radius atom compt (eg C) but dot was not masked distance=1.4 for example
                                #then, the dot is closer to nucleus of H thus get its color
                                if ( getRayon[compt2] < getRayon[compt] ):
                                    goodDots= math.sqrt((getDots[comptDots,0] - getx[compt2])**2 + (getDots[comptDots,1] - gety[compt2])**2 + (getDots[comptDots,2] - getz[compt2])**2)
                                    if ( goodDots < closestR ):
                                        closestR=goodDots
                                        mi=compt2
                                compt2=compt2+1
                        #X dummy atoms -> hole in point cloud
                        if(getA[mi]!='X'):
                            rgbi=getColor[mi]
                            getrgb[comptDots,:]=[rgb[rgbi,0], rgb[rgbi,1], rgb[rgbi,2]]
                            labeltot.append(np.vstack([getDots[comptDots], getrgb[comptDots]]))
                            if (rgbi == 1):
                                label1.append(np.vstack([getDots[comptDots], getrgb[comptDots]]))
                            elif (rgbi == 2):
                                label2.append(np.vstack([getDots[comptDots], getrgb[comptDots]]))
                            elif (rgbi == 3):
                                label3.append(np.vstack([getDots[comptDots], getrgb[comptDots]]))
                            elif (rgbi == 4):
                                label4.append(np.vstack([getDots[comptDots], getrgb[comptDots]]))
                            else:
                                print("no label for dot no %5s ?\n" %(comptDots))
                        comptDots=comptDots+1

            compt=compt+1
        
#used by both output:
#Create 10 numpy files (coordinates and color) for sensaas.py and Open3D
    if(verbose==1):
        dotsFic=open('dots.xyzrgb', 'w')
        dotslabel1=open('dotslabel1.xyzrgb', 'w')
        dotslabel2=open('dotslabel2.xyzrgb', 'w')
        dotslabel3=open('dotslabel3.xyzrgb', 'w')
        dotslabel4=open('dotslabel4.xyzrgb', 'w')
        dotspdb=open('dots.pdb', 'w')

    #print("all= %s label1= %s label2= %s label3= %s label4= %s" % (len(labeltot),len(label1),len(label2),len(label3),len(label4)))
    getDots=np.empty(shape=[len(labeltot),3], dtype='float64')
    getrgb=np.empty(shape=[len(labeltot),3], dtype='float64')
    getDots1=np.empty(shape=[len(label1),3], dtype='float64')
    getrgb1=np.empty(shape=[len(label1),3], dtype='float64')
    getDots2=np.empty(shape=[len(label2),3], dtype='float64')
    getrgb2=np.empty(shape=[len(label2),3], dtype='float64')
    getDots3=np.empty(shape=[len(label3),3], dtype='float64')
    getrgb3=np.empty(shape=[len(label3),3], dtype='float64')
    getDots4=np.empty(shape=[len(label4),3], dtype='float64')
    getrgb4=np.empty(shape=[len(label4),3], dtype='float64')

    debnom="mol"
    occup=1
    bfactor=1
    compt=0
    while(compt < len(labeltot)):
        getDots[compt]= labeltot[compt][0]
        getrgb[compt]= labeltot[compt][1]
        if(verbose==1):
            tx='%4.2f' % getDots[compt,0]
            ty='%4.2f' % getDots[compt,1]
            tz='%4.2f' % getDots[compt,2]
            dotsFic.write('%8s'%tx+'%8s'%ty+'%8s'%tz+espace5+'%5s'%getrgb[compt,0]+'%5s'%getrgb[compt,1]+'%5s'%getrgb[compt,2]+'\n')
            if(getrgb[compt,0] > 0.5 and getrgb[compt,1] > 0.5 and getrgb[compt,2] > 0.5):
                elt="H"
            elif(getrgb[compt,0] > 0.5 and getrgb[compt,1] < 0.5 and getrgb[compt,2] < 0.5):
                elt="O"
            elif(getrgb[compt,0] < 0.5 and getrgb[compt,1] > 0.5 and getrgb[compt,2] < 0.5):
                elt="C"
            else:
                elt="N"
            dotspdb.write("HETATM%5s %4s %3s     1    %8s%8s%8s  %4s %5s\n" %((compt+1),elt,debnom,tx,ty,tz,occup,bfactor))
        compt=compt+1

    compt=0
    while(compt < len(label1)):
        getDots1[compt]= label1[compt][0]
        getrgb1[compt]= label1[compt][1]
        if(verbose==1):
            tx='%4.2f' % getDots1[compt,0]
            ty='%4.2f' % getDots1[compt,1]
            tz='%4.2f' % getDots1[compt,2]
            dotslabel1.write('%8s'%tx+'%8s'%ty+'%8s'%tz+espace5+'%5s'%getrgb1[compt,0]+'%5s'%getrgb1[compt,1]+'%5s'%getrgb1[compt,2]+'\n')
        compt=compt+1

    compt=0
    while(compt < len(getDots2)):
        getDots2[compt]= label2[compt][0]
        getrgb2[compt]= label2[compt][1]
        if(verbose==1):
            tx='%4.2f' % getDots2[compt,0]
            ty='%4.2f' % getDots2[compt,1]
            tz='%4.2f' % getDots2[compt,2]
            dotslabel2.write('%8s'%tx+'%8s'%ty+'%8s'%tz+espace5+'%5s'%getrgb2[compt,0]+'%5s'%getrgb2[compt,1]+'%5s'%getrgb2[compt,2]+'\n')
        compt=compt+1

    compt=0
    while(compt < len(getDots3)):
        getDots3[compt]= label3[compt][0]
        getrgb3[compt]= label3[compt][1]
        if(verbose==1):
            tx='%4.2f' % getDots3[compt,0]
            ty='%4.2f' % getDots3[compt,1]
            tz='%4.2f' % getDots3[compt,2]
            dotslabel3.write('%8s'%tx+'%8s'%ty+'%8s'%tz+espace5+'%5s'%getrgb3[compt,0]+'%5s'%getrgb3[compt,1]+'%5s'%getrgb3[compt,2]+'\n')
        compt=compt+1

    compt=0
    while(compt < len(getDots4)):
        getDots4[compt]= label4[compt][0]
        getrgb4[compt]= label4[compt][1]
        if(verbose==1):
            tx='%4.2f' % getDots4[compt,0]
            ty='%4.2f' % getDots4[compt,1]
            tz='%4.2f' % getDots4[compt,2]
            dotslabel4.write('%8s'%tx+'%8s'%ty+'%8s'%tz+espace5+'%5s'%getrgb4[compt,0]+'%5s'%getrgb4[compt,1]+'%5s'%getrgb4[compt,2]+'\n')
        compt=compt+1
    
    if(verbose==1):
        dotsFic.close()
        dotslabel1.close()
        dotslabel2.close()
        dotslabel3.close()
        dotslabel4.close()
        dotspdb.close()

    return getDots, getrgb, getDots1, getrgb1, getDots2, getrgb2, getDots3, getrgb3, getDots4, getrgb4

############################################
