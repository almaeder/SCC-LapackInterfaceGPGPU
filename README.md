## LapackInterface
Matrix classes with internal storage that facilitate the direct invocation of LAPACK routines from C++. 

The matrix classes are primarily container classes without high level functionality incorporated as member functions. High level functionality is supplied by classes that utilized specific LAPACK routines.

There is no associated vector class, as the classes are meant to be utility classes for invoking LAPACK routines. However, since LapackMatrix and LapackMatrixCmplx16 have the standard algebraic operator member functions and provide single index element access for column matrices, instances of these classes can certainly be used as a vector class if needed. 
 
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









