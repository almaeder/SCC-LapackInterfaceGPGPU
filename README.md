## LapackInterface
Matrix classes with internal storage that facilitate the direct invocation of LAPACK routines from C++. 

The matrix classes are primarily container classes without high level functionality incorporated as member functions. High level functionality is supplied by separate classes that call specific LAPACK routines.

There is no associated vector class as the classes are meant to be utility classes for invoking LAPACK routines. However, since LapackMatrix and LapackMatrixCmplx16 have the standard algebraic operator member functions and provide single index element access for column matrices, instances of these classes declared as N X 1 matrices can be used as a vector class.  
 
### Prerequisites
C++11
### Versioning
Release : 2.0.0
### Date
Dec. 9, 2023
### Authors
Chris Anderson
### License
Lesser GPLv3  For a copy of the GNU General Public License see <http://www.gnu.org/licenses/>.
### Acknowledgements









