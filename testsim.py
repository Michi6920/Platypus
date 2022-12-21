from netgen.occ import *
from ngsolve import *
import netgen.gui
import math as m
from os import path
import numpy as np
import time
import Rotate as ro



#Test 


def calculate_torque(R_in):
    mu0 = 4e-7*pi #[kgm/A^2s^2] Permeablity of Vacuum
    
    VTK=False
    Write=False

    #Permanent Magnet
    Br=1.34 #[T]
    l_m=30e-3; b_m=10e-3;D_m=5e-3
 

    weight_Mag=11.25*1e-3 #[kg]

    #coil
    curr_dens=10*1e6 #A/m^2
    D_c=8e-3
  
    down_abs=3e-3;up_abs=3e-3
  

    #u_density=8.96 

    #Global Parameters
    Hcurl_ord=2
    H1_ord=2


   


    #Air Gap 
    A_gap=0.5e-3
   


    #Sizebox=200e-3

    Sizebox=500*1e-3



    #Stabilization Term
    Stab=1e-3


    HMAX_MAG=0.1
    HMAX_COIL=0.1
  
    if Write:
        file_path_Torque = path.relpath("Results\Sim.15.12.22_Trq_Method_Johannes\TORQUE\ full_sim_torque_DATA.txt")
        file_path_Force = path.relpath("Results\Sim.15.12.22_Trq_Method_Johannes\FORCE\ full_sim_force_DATA.txt")
        if True:
            with open(file_path_Torque, 'a') as file_TORQUE:
                file_TORQUE.write("R_in[0] , NPp ,red_Numb , mesh.ne , geo_time, mesh_time, solve_time , r_rot , Sizebox , w_c , D_c , D_m , hmax_mag, HMAX_COIL, A_gap , stab_term , Tx , Ty , Tz\n")
                with open(file_path_Force, 'a') as file_FORCE:
                    file_FORCE.write("R_in[0] , NPp ,red_Numb , mesh.ne , geo_time, mesh_time, solve_time , r_rot , Sizebox , w_c ,  D_c , D_m , hmax_mag, HMAX_COIL, A_gap , stab_term , Fx , Fy , Fz\n")
    chamfer=0.0002
    d=1e-3
    alpha=m.atan((b_m/2)/R_in[0])
    r_out=m.sqrt(R_in[0]**2+(b_m/2)**2)
    NPp=m.trunc(m.pi/(m.asin(d/m.sqrt(R_in[0]**2+(b_m/2)**2))+2*m.atan((b_m/2)/R_in[0])))
    d_gap=m.sin((2*m.pi)/(2*NPp)-2*m.atan((b_m/2)/R_in[0]))*m.sqrt(R_in[0]**2+(b_m/2)**2)
    gamma=m.asin(d_gap/R_in[0])
    beta=2*alpha+gamma

    #Calculation of max coil width
    W_c=5e-3
    if False:
        w_c=int(m.trunc(sin(beta/2)*(R_in[0]*1e3-down_abs*1e3)))*2e-3 #maximum constraint
        w_c=w_c-1e-3-d_gap
        W_c=np.arange(4,m.trunc(w_c*1e+3),1)*1e-3

    #Percentage of Rotation of Magnets
    #R_rot=[-0,-0.1,-0.2,-0.3,-0.4,-0.5]
    r_rot=0
    if False:
        R_rot=(1/4,1/2)
        R_rot=np.array(R_rot)



    print("NPp:", NPp)
    print("d_gap:", d_gap)
    print("Coil width:", W_c)


    
    
    with TaskManager():
        #Functions
        
        Force_Veks=[]
        Torque_Veks=[]
        r_Vek=np.array([0,0,R_in[0]+l_m/2])
        r_Veks=[]
        
        #Coil
        coil1_Vek=np.array([0,0,curr_dens])
        coil3_Vek=ro.rotate(-beta*180 / m.pi,coil1_Vek,)*(-1)

        if True:

            coil1_lay_2_Vek=ro.rotate(-beta*180 / m.pi*1/2,coil1_Vek)

            coil3_lay_2_Vek=ro.rotate(-beta*180 / m.pi,coil1_lay_2_Vek)*(-1)


        #Mag
        M_Vek1=CoefficientFunction((-Br/(mu0*1.05),0,0))
        M_Vek2=CoefficientFunction((Br/(mu0*1.05),0,0))
        M_Veks=[]



        #Side Coils
        coil1=Box(Pnt(D_m/2+A_gap,-W_c/2,R_in[0]-down_abs), Pnt(D_m/2+A_gap+D_c,W_c/2,R_in[0]+l_m+up_abs))
        coil1 =coil1.MakeChamfer (coil1.edges,chamfer)
        coil3=coil1.Rotate(Axis((0,0,0), X),-beta*180 / m.pi)

        #Second Layer Coils
        coil1_lay_2=coil1.Mirror(Axis((0,0,0), Z))
        
        coil1_lay_2=coil1_lay_2.Rotate(Axis((0,0,0), X),-beta*180 / m.pi*1/2)
        
            
        coil3_lay_2=coil1_lay_2.Rotate(Axis((0,0,0), X),-beta*180 / m.pi)
    

        #Magnet
        mag0 =Box(Pnt(-D_m/2,-b_m/2,R_in[0]), Pnt(D_m/2,b_m/2,R_in[0]+l_m))
        mag0 = mag0.MakeChamfer (mag0.edges,chamfer)
        mag0= mag0.Rotate(Axis((0,0,0), X),r_rot*beta/m.pi*180)

        #Iterations
        rot_mag=-(beta)* 180 / m.pi*np.arange(0,2*NPp,1)
        rot_coil=-(2*beta)* 180 / m.pi*np.arange(0,NPp,1)


        #Air Box
        box = Box(Pnt(-Sizebox,-Sizebox,-Sizebox),Pnt(Sizebox,Sizebox,Sizebox))

        #Setup Size
        red_array=4
        if False:
            if NPp%2==1:
                red_array=[2*NPp]
            else:
                red_array=[2*NPp,NPp,int(NPp/2)]
                i=0
                while i < len(red_array):
                    if red_array[i]<2*NPp:
                        if red_array[i]%2 ==1:
                            red_array[i]+=3
                        else:
                            red_array[i]+=2
                    i+=1

        
        start_time_total=time.time()
        start_time_geo=time.time()
        
        
            
        if Write:
            file_path_DATA = path.relpath("Results\Sim.15.12.22_Trq_Method_Johannes\DATA\Curl2 H2 {0}, {1} , {2}, {3} , {4} , {5} , {6} ,{7} , {8} , {9} , {10} , {11} .txt".format(R_in[0], NPp ,red_array[b] , r_rot , Sizebox , w_c , D_c ,  D_m , hmax_mag, HMAX_COIL, A_gap , stab))
            file_path_INFO =path.relpath("Results\Sim.15.12.22_Trq_Method_Johannes\INFO\Curl2 H2 {0}, {1} , {2} , {3} , {4} , {5} , {6} ,{7} , {8} , {9} , {10} , {11} .txt".format(R_in[0], NPp ,red_array[b] , r_rot , Sizebox , w_c , D_c ,  D_m , hmax_mag, HMAX_COIL, A_gap , stab))
        if Write:     
            with open(file_path_DATA, 'a') as file_DATA:
                file_DATA.write("Nu-Mag: %s ,R_in[0], NPp ,red_Numb , r_rot , Sizebox , w_c ,  D_c ,  D_m ,  hmax_mag, HMAX_COIL, A_gap , stab ,Fx,Fy,Fz,Tx,Ty,Tz\n" %(red_array[b])) 
        red_Numb=red_array
        
        rot1_mag=rot_mag[0:red_Numb]
        rot1_coil=rot_coil[0:int(red_Numb/2)]
        
        for u in range(len(rot1_coil)):
            M_Veks.append(M_Vek1)
            M_Veks.append(M_Vek2)
        bodies=[]
        coil1s=[]
        
        coil3s=[]
        
        coil1_lay_2s=[]
        coil3_lay_2s=[]
        
        

        coil1_Veks=[]
        coil3_Veks=[]
        coil1_lay_2_Veks=[]
        coil3_lay_2_Veks=[]

        i=0
        for r in rot1_mag:

            body_rot=mag0.Rotate(Axis((0,0,0), X),r).bc('gamma'+str(i))
            body_rot.maxh=HMAX_MAG
            body_rot.mat("Mag"+str(i))
            r_Vek_rot=ro.rotate(r,r_Vek)
            r_Veks.append(r_Vek_rot)
        
            if (i%2)==0:
                coil1_rot=coil1.Rotate(Axis((0,0,0), X),r)
                
                coil3_rot=coil3.Rotate(Axis((0,0,0), X),r)

                coil1_lay_2_rot=coil1_lay_2.Rotate(Axis((0,0,0), X),r)

                coil3_lay_2_rot=coil3_lay_2.Rotate(Axis((0,0,0), X),r)
                
                coil1_rot.maxh=HMAX_COIL

                coil3_rot.maxh=HMAX_COIL
                
                coil1_lay_2_rot.maxh=HMAX_COIL

                coil3_lay_2_rot.maxh=HMAX_COIL
                




                coil1_rot.mat(str(i)+"coil1_vec")
            
                coil3_rot.mat(str(i)+"coil3_vec")

                coil1_lay_2_rot.mat(str(i)+"coil1_lay2_vec")
            
                coil3_lay_2_rot.mat(str(i)+"coil3_lay2_vec")
                

                #alternating Vectors
                coil1_Veks.append(ro.rotate(r,coil1_Vek))
                coil3_Veks.append(ro.rotate(r,coil3_Vek))
                coil1s.append(coil1_rot)
                coil3s.append(coil3_rot)

                coil1_lay_2_Veks.append(ro.rotate(r,coil1_lay_2_Vek))
                coil3_lay_2_Veks.append(ro.rotate(r,coil3_lay_2_Vek))
                coil1_lay_2s.append(coil1_lay_2_rot)
                coil3_lay_2s.append(coil3_lay_2_rot)

            bodies.append(body_rot)
            
    
            i=i+1
    

        
        for c in range(len(coil1s)):
            bodies.append(coil1s[c])
            bodies.append(coil3s[c])

            
            bodies.append(coil1_lay_2s[c])
            bodies.append(coil3_lay_2s[c])
            if False:
                bodies.append(coil4s[c])
                bodies.append(coil2s[c])




        y=0
        air=box
    
        while y<len(bodies):
            
            air=air-bodies[y]
            
            y=y+1
        air.mat("air")
        bodies.append(air)
    
        mur = {"air" : 1.00000037}
        for co in range(len(rot1_mag)):
            x="Mag"+str(co)
            if co%2==0:
                x1=str(co)+"coil1_vec"
                x3=str(co)+"coil3_vec"
                x1_2lay=str(co)+"coil1_lay2_vec"
                x3_2lay=str(co)+"coil3_lay2_vec"
                mur[x1]=0.999994
                mur[x3]=0.999994
                mur[x1_2lay]=0.999994
                mur[x3_2lay]=0.999994
    
            mur[x]=1.05
            
        
    
    
        geo = OCCGeometry(bodies)  
        if True:
            geo_time=(time.time()-start_time_geo)/60  
            start_time_mesh=time.time()
            mesh = Mesh(geo.GenerateMesh(maxh=1))
        
            mesh_time=(time.time()-start_time_mesh)/60 
            print(mesh.ne)
            mesh_temp=mesh.ne
            if False:
                if(mesh.ne> 1410000):
                    print("Assembly will fail")
                    if Write:
                        with open(file_path_INFO, 'a') as file_INFO: 
                            file_INFO.write ("Numb_Mags:%s\n"%(red_array[b]))
                            file_INFO.write("Mesh_Elements:%s\n"%(mesh.ne))
                            file_INFO.write ("Assembly will fail\n")
                            file_INFO.write ("R_in[0]:%s coil width:%s\n"%(R_in[0],w_c))
                    b+=1
                    continue
        
            start_time_solve=time.time()
            if Write:
                with open(file_path_INFO, 'x') as file_INFO:
                    file_INFO.write ("A_gap:%s w_c:%s   D_c:%s\n"%(A_gap,w_c,D_c))
            
            if True:
                if True:
                    fes = HCurl(mesh, order=Hcurl_ord, nograds=True)
                    u,v = fes.TnT()
                    nu = 1.0 / mu0 / CoefficientFunction( [mur[mat] for mat in mesh.GetMaterials()] )
                    
                    


                    a = BilinearForm(fes)
                    a += SymbolicBFI (nu*curl(u)*curl(v) + nu*Stab*u*v)
                    c = Preconditioner(a, "bddc")

                    f = LinearForm(fes)
                    
                    
                    if True:
                        for co in range(len(rot1_mag)):

                            f += SymbolicLFI(M_Veks[co] * curl(v), definedon=mesh.Materials("Mag"+str(co)))
                        
                        for co in range(len(rot1_coil)):

                            f += SymbolicLFI(CoefficientFunction((coil1_Veks[co][0],coil1_Veks[co][1],coil1_Veks[co][2]))*v, definedon=mesh.Materials((str(co*2)+'coil1_vec')))


                            f += SymbolicLFI(CoefficientFunction((coil3_Veks[co][0],coil3_Veks[co][1],coil3_Veks[co][2]))*v, definedon=mesh.Materials((str(co*2)+"coil3_vec")))


                            f += SymbolicLFI(CoefficientFunction((coil1_lay_2_Veks[co][0],coil1_lay_2_Veks[co][1],coil1_lay_2_Veks[co][2]))*v, definedon=mesh.Materials((str(co*2)+"coil1_lay2_vec")))


                            f += SymbolicLFI(CoefficientFunction((coil3_lay_2_Veks[co][0],coil3_lay_2_Veks[co][1],coil3_lay_2_Veks[co][2]))*v, definedon=mesh.Materials((str(co*2)+"coil3_lay2_vec")))


                            if False:
                                f += SymbolicLFI(curr_dens*gtau*v, definedon=mesh.Materials(("coil4_vec")))
                                f += SymbolicLFI(curr_dens*ftau*v, definedon=mesh.Materials(("coil2_vec")))
                        
                    
                    
                    
                    
                    a.Assemble()
                    f.Assemble()
                    gfu = GridFunction(fes)
                    Kinv = CGSolver(a.mat, c.mat)
                    gfu.vec.data = Kinv * f.vec
                    Draw(curl(gfu), mesh, "B-field")
                    solve_time=(time.time()-start_time_solve)/60 
                    start_time_vtk=time.time()
                    if VTK:
                        
                        vtk = VTKOutput(ma=mesh,coefs=[curl(gfu)],names=["B-Field"],filename="vtk_Exceptions_"+str(red_array[b])+str(R_in[0])+str(w_c*1e3),subdivision=0,legacy=True)
                        vtk.Do()
                        
                    vtk_time=(time.time()-start_time_vtk)/60 
                    start_time_MST=time.time()
                    if True:
                    # Maxwell stress tensor ...
                        B = curl(gfu)
                        sigma = nu * (OuterProduct(B, B) - 0.5 * (B*B) * Id(3))

                        # compute force in H1-space:
                        # compute f = div sigma  in weak form
                        fes = VectorH1(mesh, order=H1_ord)
                        u,v = fes.TnT()
                        force = LinearForm(fes)
                        force += SymbolicLFI( InnerProduct(sigma, grad(v)))
                        force.Assemble()

                        mass = BilinearForm(fes)
                        mass += SymbolicBFI(u*v, BND)
                        mass.Assemble()
                        MST_time=(time.time()-start_time_MST)/60 
                        Fx=0
                        start_time_Force_torque=time.time()
                        Force_tot=np.array([0,0,0])
                        Torque_tot=np.array([0,0,0])
                        Force_tot_red=np.array([0,0,0])
                        Torque_tot_red=np.array([0,0,0])

                        
                        for co in range(len(rot1_mag)):

                            gfforce = GridFunction(fes)
                            gfforce.vec.data = mass.mat.Inverse(inverse="sparsecholesky", freedofs=fes.GetDofs(mesh.Boundaries('gamma'+str(co)))) * force.vec
                            F_mst =Integrate(gfforce, mesh, BND, definedon=mesh.Boundaries('gamma'+str(co)))
                            T_mst=np.cross(r_Veks[co],F_mst)
                            if co>0 and co<len(rot1_mag)-1:
                                Force_tot_red=Force_tot_red+F_mst
                                Torque_tot_red=Torque_tot_red+T_mst
                            Force_tot=Force_tot+F_mst
                            Torque_tot=Torque_tot+T_mst

                        Force_tot_red=Force_tot_red*3.6
                        Torque_tot_red=Torque_tot_red*3.6
    return (Torque_tot_red[0]-1)

                                            
                                                    
                                                                
                                                                                
                    

