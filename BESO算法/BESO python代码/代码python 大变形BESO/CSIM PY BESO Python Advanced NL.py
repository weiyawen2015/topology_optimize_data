"""General 3D topology optimization code by Zhihao Zuo. Note that the CAE file shall contain a model 'Model-1' with a dependent part 'Part-1' and a static step 'Step-1'."""
import math,customKernel
from abaqus import getInput,getInputs
from odbAccess import openOdb
## Function of formatting Abaqus model for stiffness optimisation (nonlinear)
def fmtMdb(Mdb):
    mdl = Mdb.models['Model-1']
    part = mdl.parts['Part-1']
    # Build sections and assign solid section
    mdl.Material('Material01').Elastic(((1.0, 0.3), ))
    mdl.HomogeneousSolidSection('sldSec','Material01')
    mdl.Material('Material02').Elastic(((0.001**3, 0.3), ))
    mdl.HomogeneousSolidSection('voidSec','Material02')
    part.SectionAssignment(part.Set('ss',part.elements),'sldSec')
    # Define output request
    mdl.FieldOutputRequest('SEDensity','Step-1',variables=('ENER', ))
    mdl.HistoryOutputRequest('ExtWork','Step-1',variables=('ALLWK', ))
## Function of running FEA for for raw sensitivities and objective function (nonlinear)
def FEA(Iter,Mdb,Xe,Ae):
    try:
        Mdb.Job('Design_Job'+str(Iter),'Model-1').submit()
        Mdb.jobs['Design_Job'+str(Iter)].waitForCompletion()
    except AbaqusException, message:
        print "Error occured: ", message
        sys.exit(1)
    opdb = openOdb('Design_Job'+str(Iter)+'.odb')
    seng = opdb.steps['Step-1'].frames[-1].fieldOutputs['SENER'].values
    for en in seng: Ae[en.elementLabel]=en.data/Xe[en.elementLabel]
    seng = opdb.steps['Step-1'].frames[-1].fieldOutputs['PENER'].values
    for en in seng: Ae[en.elementLabel]+=en.data/Xe[en.elementLabel]
    obj=opdb.steps['Step-1'].historyRegions['Assembly ASSEMBLY'].historyOutputs['ALLWK'].data[-1][1]
    opdb.close()
    return obj
## Function of preparing filter map (Fm = {elm1:[[el1,el2,...],[wf1,wf2,...]],...}) (NumPy)
def preFlt(Rmin,Elmts,Nds,Fm):
    import numpy as np
    # Calculate element centre coordinates
    elm, c0 = np.zeros(len(Elmts)), np.zeros((len(Elmts),3))
    for i in range(len(elm)):
        elm[i] = Elmts[i].label
        nds = Elmts[i].connectivity
        for nd in nds: c0[i] = np.add(c0[i],np.divide(Nds[nd].coordinates,len(nds)))
    # Weighting factors
    for i in range(len(elm)):
        Fm[elm[i]] = [[],[]]
        for j in range(len(elm)):
            dis = np.square(np.sum(np.power(np.subtract(c0[i],c0[j]),2)))
            if dis<Rmin: 
                Fm[elm[i]][0].append(elm[j])
                Fm[elm[i]][1].append(Rmin - dis)
        Fm[elm[i]][1] = np.divide(Fm[elm[i]][1],np.sum(Fm[elm[i]][1]))
## Function of filtering sensitivities
def fltAe(Ae,Fm):
    raw = Ae.copy()
    for el in Fm.keys():
        Ae[el] = 0.0
        for i in range(len(Fm[el][0])): Ae[el]+=raw[Fm[el][0][i]]*Fm[el][1][i]
## Function of optimality update for design variables and Abaqus model
def BESO(Vf,Xe,Ae,Part,Elmts):
    lo, hi = min(Ae.values()), max(Ae.values())
    tv = Vf*len(Elmts)
    while (hi-lo)/hi > 1.0e-5:
        th = (lo+hi)/2.0
        for key in Xe.keys(): Xe[key] = 1.0 if Ae[key]>th else 0.001
        if sum(Xe.values())-tv>0: lo = th
        else: hi = th
    # Label elements as solid or void
    vlb, slb = [], []
    for el in Elmts:
        if Xe[el.label] == 1.0: slb.append(el.label)
        else: vlb.append(el.label)
    # Assign solid and void elements to each section
    Part.SectionAssignment(Part.SetFromElementLabels('ss',slb),'sldSec')
    Part.SectionAssignment(Part.SetFromElementLabels('vs',vlb),'voidSec')
## ====== MAIN PROGRAM ======
if __name__ == '__main__':
    # Set parameters and inputs
    pars = (('VolFrac:','0.5'), ('Rmin:', '1'), ('ER:', '0.02'))
    vf,rmin,ert = [float(k) if k!=None else 0 for k in getInputs(pars,dialogTitle='Parameters')]
    if vf<=0 or rmin<0 or ert<=0: sys.exit()
    mddb = openMdb(getInput('Input CAE file:',default='Test.cae'))
    # Design initialization
    fmtMdb(mddb)
    part = mddb.models['Model-1'].parts['Part-1']
    elmts, nds = part.elements, part.nodes
    oh, vh = [], []
    xe, ae, oae, fm = {}, {}, {}, {}
    for el in elmts: xe[el.label] = 1.0
    if rmin>0: preFlt(rmin,elmts,nds,fm)
    # Optimisation iteration
    change, iter, obj = 1, -1, 0
    while change > 0.001:
        iter += 1
        # Run FEA
        oh.append(FEA(iter,mddb,xe,ae))
        # Process sensitivities
        if rmin>0: fltAe(ae,fm)
        if iter > 0: ae=dict([(k,(ae[k]+oae[k])/2.0) for k in ae.keys()])
        oae = ae.copy()
        # BESO optimisation
        vh.append(sum(xe.values())/len(xe))
        nv = max(vf,vh[-1]*(1.0-ert))
        BESO(nv,xe,ae,part,elmts)
        if iter>10: change=math.fabs((sum(oh[iter-4:iter+1])-sum(oh[iter-9:iter-4]))/sum(oh[iter-9:iter-4]))
    # Save results
    mddb.customData.History = {'vol':vh,'obj':oh}
    mddb.saveAs('Final_design.cae')
