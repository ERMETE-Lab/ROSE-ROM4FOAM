import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Default values
DeltaT = 1.
g = 9.81

parallel = True
num_proc = 4

# Get Current Directory
w_dir = os.getcwd()

# Train/Test Parameters
dRe = 5.
dRi = 0.4

train = False
if train:
    print('Creating Train Snapshots')
    Re = np.arange(15,  150+dRe/2, dRe) # [15, 150]
    Ri = np.arange(0.2,   5+dRi/2, dRi) # [0.2, 1., 5.]
    if not os.path.exists(w_dir+'/TrainSet'):
        os.system("mkdir "+w_dir+'/TrainSet')
else:
    print('Creating Test Snapshots')
    Re = np.arange(15+dRe/2,  150+dRe/2, dRe) 
    Ri = np.arange(0.2+dRi/2,   5+dRi/2, dRi)
    # Re = [15, 150]
    # Ri = [0.2, 1., 5.]
    if not os.path.exists(w_dir+'/TestSet'):
        os.system("mkdir "+w_dir+'/TestSet')

def clean_case():
    os.system("rm -r Case_*")

caseI = 0
parameters = np.zeros((len(Re)*len(Ri), 2))


fig_path = w_dir+"/Residual_Fig"
if not os.path.exists(fig_path):
    os.system("mkdir "+fig_path)

for ReI in Re:
    for RiI in Ri:

        print('Case solving for Re = {:.2f}'.format(ReI)+' and Ri = {:.2f}'.format(RiI))

        # Copy base case to new folder
        if train:
            folder_name = w_dir+"/TrainSet/"+f'Case_{caseI+0:03}_Re{ReI:.2f}_Ri{RiI:.2f}'
        else:
            folder_name = w_dir+"/TestSet/"+f'Case_{caseI+0:03}_Re{ReI:.2f}_Ri{RiI:.2f}'
        os.system("cp -r "+w_dir+"/BaseCase "+folder_name)

        # Set path
        path = folder_name
        os.chdir(path)

        # Update Readme with new Re and Ri
        with open(path+"/README.md", "r") as f:
            lines = f.readlines()
            lines[0] = 'Case solved for Re = {:.2f}'.format(ReI)+' and Ri = {:.2f}'.format(RiI)
        with open(path+"/README.md", "w") as f:
            f.writelines(lines)

        # Read constant/transportProperties file
        with open(path+"/constant/transportProperties", "r") as f:
            lines = f.readlines()

        # Change lines for nu and beta (i.e. Re and Ri)
        lines[21] = 'nu              [0 2 -1 0 0 0 0] {:.8f}'.format(1/ReI)+';\n'
        lines[24] = 'beta            [0 0 0 -1 0 0 0] {:.8f}'.format(RiI/g/DeltaT)+';\n'

        # Write updated constant/transportProperties file
        with open(path+"/constant/transportProperties", "w") as f:
            f.writelines(lines)

        # Generate the mesh
        os.system("blockMesh > log.blockMesh")
        os.system("renumberMesh -overwrite > log.renumber")

        # Run the solver either in parallel or in series
        if parallel:
            os.system("decomposePar > log.decompose")
            # os.system("mpirun -np "+str(num_proc)+" buoyantBoussinesqSimpleFoam")
            os.system("mpirun -np "+str(num_proc)+" buoyantBoussinesqSimpleFoam > log.solve" )
            os.system("reconstructPar > log.reconstruct")
            os.system("rm -r proc*")
        else:
            os.system("buoyantBoussinesqSimpleFoam > log.solve")

        # Plot Residuals
        data = pd.read_csv(path+"/postProcessing/residuals/0/residuals.dat",skiprows=1, delimiter='\s+').iloc[:, 1:].shift(+1,axis=1).drop(["Time"], axis= 1)
        plot = data.plot(logy= True, logx=False, figsize=(6,4))
        fig = plot.get_figure()
        ax = plt.gca()
        ax.legend(loc='upper right')
        ax.set_xlabel("Iterations")
        ax.set_ylabel("Residuals")
        plt.grid()
        plt.tight_layout()
        plt.savefig(fig_path+f"/residuals_case_{caseI+0:03}.pdf", format='pdf', dpi=300, bbox_inches='tight')
        plt.close()
        parameters[caseI, :] = [ReI, RiI]
        caseI += 1
        
np.savetxt(w_dir+'/parameters.txt', parameters, delimiter=',')

folder_list = list()
caseI = 0
for ReI, RiI in parameters[:, :]:

        if train is False:
            folder_name =f'"TestSet/Case_{caseI+0:03}_Re{ReI:.2f}_Ri{RiI:.2f}"\n'
        else:
            folder_name =f'"Case_{caseI+0:03}_Re{ReI:.2f}_Ri{RiI:.2f}"\n'
        folder_list.append(folder_name)
        caseI += 1

if train:
    with open(w_dir+'/TrainSet/train_folders.txt', "w") as f:
        f.writelines(folder_list)
else:
    with open(w_dir+'/TestSet/test_folders.txt', "w") as f:
        f.writelines(folder_list)