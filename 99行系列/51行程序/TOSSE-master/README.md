# TOSSE
An Efficient 51 lines Matlab code for topology optimization.

TOSSE (Topology Optimization Same Size Elements) is a Matlab code for 2D and 3D topology design problems. 
The code uses the classic 88 lines code known as TOP88 as basis in order to develop an hard 0-1 evolutionary algorithm that, at every iteration, hard kills the elements. The new code is comprised of 51 lines without sacrificing any readability, making it useful for practitioners who want to approach the field.

The algorithm shows an efficiency on averange superior to TOP88 and structures with almost none checkered board patterns.

For more details on the theory and the numerical results you can check the paper available on:
https://arxiv.org/abs/1902.00877

## Usage
On this project three codes are available: 
```
tosse.m
tosse_cant.m
tosse3d.m
```
The first is a topology optimization code for the Messerschmitt-Bolkow-Blohm (MBB) beam. The code can be launched by typing in the Matlab terminal the following command:
```
tosse(nelx,nely,volfrac,mu)
```
where ```nelx``` is the number of elements on the x-axis, ```nely``` is the number of elements on the y axis, ```volfrac``` is the desired volume in the final design and ```mu``` is the volume reduction parameter.

A practical example of a call is:
```
tosse(180,60,0.5,0.97)
```
For a 180X60 design domain with half of the original volume.
### tosse_cant.m
This code solves 2D topologies for the cantilever beam. The code can be launched in a similar fashion to tosse.m:
```
tosse_cant(nelx,nely,volfrac,mu,sym)
```
where the parameter ```sym``` indicates if the algorithm has to enforce symmetry in the final design.
### tosse3d.m
This code solved 3D topologies for the cantilever beam. The launch command on the terminal is:
```
tosse3d(nelx,nely,nelz,volfrac,vc,mu)
```
where ```nelz``` is the number of elements on the z-axis.

## Extensions
This code can be considered as a basis for future more complex codes in topology optimization. Many of the extensions that can be found in literature like other boundary conditions, passive elements, displacement constraints and so on (see for example [1](https://link.springer.com/article/10.1007/s00158-010-0594-7),[2](https://link.springer.com/article/10.1007/s00158-010-0487-9),[3](https://link.springer.com/article/10.1007/s00158-014-1107-x)) can be applied with small modifications to this code as well.
## License
The software is distribuited under the GNU License. Check the License.md file for more information.
