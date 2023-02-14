/* Eric Joyce, Stevens Institute of Technology, 2019; ejoyce@stevens.edu

   19dec19:  Subsamples the given mesh as a point-cloud.
             Assumes vertices of the given mesh have attributes x, y, z, r, g, b, and isCam,
             but does NOT include isCam in the constructed point-cloud.
*/

#include <fstream>
#include <iostream>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <Eigen/Dense>

/*
#define __PSAMPLE_DEBUG 1
*/
/*
#define __PSAMPLE_LOAD_DEBUG 1
*/

#define SAMPLE_MAIN_MALLOC_MAX 16777215                             /* Largest number of samples to accumulate before caching them */

typedef struct VertexType
  {
    float x;                                                        //  The specifications of these data types
    float y;                                                        //  come from the ETH PLY file headers.
    float z;                                                        //  Even if we do not use all of these for
    unsigned char r;                                                //  our purposes, we need to know how large
    unsigned char g;                                                //  they are so that we can skip past them
    unsigned char b;                                                //  in the binary data.
    bool cam;                                                       //  uchar in PLY spec: convert to bool
  } Vertex;

typedef struct FaceType
  {
                                                                    //  An initial value, n, is always going to be = 3. Ignore it.
    unsigned int v[3];                                              //  Indices of Vertices making up this Face
    double barycenter[3];                                           //  X, Y, and Z components of the Face barycenter
    double normal[3];                                               //  X, Y, and Z components of unit normal to the Face
    double area;                                                    //  Area of this face
  } Face;

/*
g++ -Wall -I ./ psample.cpp -lm -o psample
./psample terrains.mincut.KB.ply 0.0001
*/

using namespace std;
using namespace Eigen;

void writePLY(const char*, Vertex*, unsigned long, unsigned int);
unsigned long sampleFaces(double, Vertex*, Face*, unsigned int, Vertex**, unsigned int*);
void cacheSamples(Vertex**, unsigned long, unsigned int*);
bool loadMesh(const char*, unsigned int*, Vertex**, Face**);
bool getMeshInfo(const char*, unsigned int**);
void computeFaceBarycenter(Vertex*, Face**, unsigned int);
void computeFaceArea(Vertex*, Face**, unsigned int);
void computeFaceNormal(Vertex*, Face**, unsigned int);
float crossProduct(float*, float*, float**);
double bary2d(Vector3f, Vector3f, double*);
double min3(double, double, double);
double max3(double, double, double);
void usage(void);

int main(int argc, char* argv[])
  {
    Vertex* V;                                                      //  Array of vertices
    Face* F;                                                        //  Array of faces

    Vertex* points;                                                 //  Array of vertices derived from faces
    unsigned long len = 0L;                                         //  Length of this array

    double sampleRate = 0.01;                                       //  Scalar by which we multiply unit vectors on each Face
    unsigned int cacheCount = 0;                                    //  Number of cache files created during sampling.
    unsigned int* fileInfo;                                         //  Array of data about the given PLY file
    unsigned char fileInfoLen;                                      //  Length of that array

    fileInfoLen = 4;                                                //  There are 4 things we need to know about the PLY:
    if(argc < 2)                                                    //  [0] Number of vertices
      {                                                             //  [1] Number of faces
        usage();                                                    //  [2] Number of tetrahedra (not needed here)
        return 0;                                                   //  [3] Offset of the first byte of PLY binary data
      }                                                             //  And make sure we have enough arguments!
    if(argc > 2)                                                    //  This argument is optional:
      sampleRate = (double)atof(argv[2]);                           //  If given a 'sampleRate', use that

                                                                    //  Allocate file info array
    if((fileInfo = (unsigned int*)malloc(fileInfoLen * sizeof(int))) == NULL)
      {
        cout << "ERROR: Unable to allocate file info array." << endl;
        return 0;
      }

    if(!getMeshInfo(argv[1], &fileInfo))                            //  Get PLY header information
      {
        cout << "ERROR: Unable to allocate file info array." << endl;
        return 0;
      }

    if(!loadMesh(argv[1], fileInfo, &V, &F))                        //  Only load vertices and faces
      {
        if(fileInfo[0] > 0)
          free(V);
        if(fileInfo[1] > 0)
          free(F);
        free(fileInfo);
        return 0;
      }
                                                                    //  Sample the faces
    len = sampleFaces(sampleRate, V, F, fileInfo[1], &points, &cacheCount);

    if(fileInfo[0] > 0)                                             //  These are no longer necessary
      free(V);
    if(fileInfo[1] > 0)
      free(F);
    free(fileInfo);

    writePLY("psample.ply", points, len, cacheCount);               //  Create the new file

    if(len > 0L)
      free(points);

    return 0;
  }

/* Subsampling is complete, and we write the point cloud to file.
   In addition to the current contents of array 'V', of length 'len', check for 'cacheCtr' cache files
   and include them in the final output.
   Note that here we do NOT include the isCam attribute. */
void writePLY(const char* filename, Vertex* V, unsigned long len, unsigned int cacheCtr)
  {
    ifstream ifh;                                                   //  Input file handle (cached file)
    ofstream fh;                                                    //  Output file handle
    char buffer[64];
    unsigned long i;
    unsigned int j;
    unsigned long vtotal = len;                                     //  Total will be *at least* = 'len'
    unsigned long filecount;
    bool breakout;                                                  //  Loop-break test
    float x, y, z;                                                  //  Holders for cached vertex values
    unsigned char r, g, b;                                          //  (Always the same = 128)

    #ifdef __PSAMPLE_DEBUG
    cout << "writePLY(" << +cacheCtr << " files cached)" << endl;
    #endif

    for(j = 0; j < cacheCtr; j++)
      {
        sprintf(buffer, "psample.%d.cache", j);
        ifh.open(buffer, ios::binary);                              //  Open the cache file (for binary reading)
        ifh.seekg(0, ios::beg);                                     //  Start at the beginning
        ifh.read(reinterpret_cast<char*>(&filecount), sizeof(long));

        cout << buffer << " contributes " << filecount << " vertices" << endl;

        vtotal += filecount;
        ifh.close();
      }

    cout << "+ " << len << " vertices in accumulator array" << endl;
    cout << "TOTAL = " << vtotal << " vertices" << endl;

    fh.open(filename);                                              //  Open the file
    fh << "ply" << endl;                                            //  Ply format header
    fh << "format binary_little_endian 1.0" << endl;                //  Declare endianness
    fh << "element vertex " << +vtotal << endl;                     //  Declare number of vertices
    fh << "property float x" << endl;                               //  Vertex properties are the same as in the original PLY
    fh << "property float y" << endl;
    fh << "property float z" << endl;
    fh << "property uchar red" << endl;
    fh << "property uchar green" << endl;
    fh << "property uchar blue" << endl;
    fh << "end_header" << endl;                                     //  Close the header

    for(j = 0; j < cacheCtr; j++)
      {
        sprintf(buffer, "psample.%d.cache", j);
        ifh.open(buffer, ios::binary);                              //  Open the cache file (for binary reading)
        ifh.seekg(0, ios::beg);                                     //  Start at the beginning
        ifh.read(reinterpret_cast<char*>(&filecount), sizeof(long));//  Read this to advance the file pointer and know when count is reached.
        i = 0L;
        breakout = false;
        while(!ifh.eof() && !breakout)
          {
            ifh.read(reinterpret_cast<char*>(&x), sizeof(float));   //  Read out of cache
            ifh.read(reinterpret_cast<char*>(&y), sizeof(float));
            ifh.read(reinterpret_cast<char*>(&z), sizeof(float));
            ifh.read(reinterpret_cast<char*>(&r), sizeof(char));
            ifh.read(reinterpret_cast<char*>(&g), sizeof(char));
            ifh.read(reinterpret_cast<char*>(&b), sizeof(char));

            fh.write(reinterpret_cast<char*>(&x), sizeof(float));   //  Write into final file
            fh.write(reinterpret_cast<char*>(&y), sizeof(float));
            fh.write(reinterpret_cast<char*>(&z), sizeof(float));
            fh.write(reinterpret_cast<char*>(&r), sizeof(char));
            fh.write(reinterpret_cast<char*>(&g), sizeof(char));
            fh.write(reinterpret_cast<char*>(&b), sizeof(char));

            if(++i == filecount)
              breakout = true;                                      //  Done with this one!
          }

        ifh.close();                                                //  Close cache file
      }

    for(i = 0L; i < len; i++)                                       //  Write every Vertex in 'V' array
      {
        fh.write(reinterpret_cast<char*>(&V[i].x), sizeof(float));
        fh.write(reinterpret_cast<char*>(&V[i].y), sizeof(float));
        fh.write(reinterpret_cast<char*>(&V[i].z), sizeof(float));
        fh.write(reinterpret_cast<char*>(&V[i].r), sizeof(char));
        fh.write(reinterpret_cast<char*>(&V[i].g), sizeof(char));
        fh.write(reinterpret_cast<char*>(&V[i].b), sizeof(char));
      }

    fh.close();                                                     //  Close the file. Done!

    return;
  }

/* Given the number of points per unit area, sample all Faces in array F
   and write points to array 'points' and return the length of 'points'.
   Return 0 if something terrible happened (or if there were simply no points.) */
unsigned long sampleFaces(double sampleRate, Vertex* V, Face* F, unsigned int lenF, Vertex** points, unsigned int* cacheCount)
  {
    unsigned long len = 0;                                          //  Length of the array to be filled

    Vector3f A, B, C;                                               //  Transformed versions of vertex-0, vertex-1, and vertex-2
    Vector3f u, v, w;                                               //  Used to compute the rotation matrix
    float mag;                                                      //  Vector magnitude holder

    Matrix3f R;                                                     //  Rotation matrix

    double minX, minY;                                              //  Bounding box around triangle
    double maxX, maxY;
    double q[3];                                                    //  Hold values for triangle test
    double p[2];

    Vertex* samples;                                                //  Array of samples for a given Face
    unsigned int sampleLen = 0;                                     //  Length of the array

    unsigned int i, j;

    for(i = 0; i < lenF; i++)                                       //  For each Face
      {
                                                                    //  (This much is helpful)
        cout << "Sampling Face " << +i << "/" << +lenF << " at " << sampleRate << endl;

        #ifdef __PSAMPLE_DEBUG
        cout << "Sampling Face " << +i << "/" << +lenF << " at " << sampleRate << endl;
        cout << "  F[" << +i << "] = " << +F[i].v[0] << " " << +F[i].v[1] << " " << +F[i].v[2] << endl;
        cout << "  area = " << F[i].area << endl;
        #endif

        A(0) = 0.0;                                                 //  Subtract vertex[0] from vertex[0]
        A(1) = 0.0;
        A(2) = 0.0;

        B(0) = (double)(V[ F[i].v[1] ].x - V[ F[i].v[0] ].x);       //  Subtract vertex[0] from vertex[1]
        B(1) = (double)(V[ F[i].v[1] ].y - V[ F[i].v[0] ].y);
        B(2) = (double)(V[ F[i].v[1] ].z - V[ F[i].v[0] ].z);
        mag = sqrt(B(0) * B(0) + B(1) * B(1) + B(2) * B(2));

        C(0) = (double)(V[ F[i].v[2] ].x - V[ F[i].v[0] ].x);       //  Subtract vertex[0] from vertex[2]
        C(1) = (double)(V[ F[i].v[2] ].y - V[ F[i].v[0] ].y);
        C(2) = (double)(V[ F[i].v[2] ].z - V[ F[i].v[0] ].z);

        u(0) = B(0) / mag;                                          //  Vector u = normalized B
        u(1) = B(1) / mag;
        u(2) = B(2) / mag;

        w(0) = u(1) * C(2) - u(2) * C(1);                           //  Vector w = u cross C
        w(1) = u(2) * C(0) - u(0) * C(2);
        w(2) = u(0) * C(1) - u(1) * C(0);
        mag = sqrt(w(0) * w(0) + w(1) * w(1) + w(2) * w(2));        //  Normalize w
        w(0) /= mag;
        w(1) /= mag;
        w(2) /= mag;

        v(0) = u(1) * w(2) - u(2) * w(1);                           //  Vector v = u cross w
        v(1) = u(2) * w(0) - u(0) * w(2);
        v(2) = u(0) * w(1) - u(1) * w(0);
                                                                    //  Rotation matrix
        R(0, 0) = u(0);  R(0, 1) = v(0);  R(0, 2) = w(0);
        R(1, 0) = u(1);  R(1, 1) = v(1);  R(1, 2) = w(1);
        R(2, 0) = u(2);  R(2, 1) = v(2);  R(2, 2) = w(2);
        //  RotMat = [ u[0] v[0] w[0] ]                             Takes (1,0,0) to u
        //           [ u[1] v[1] w[1] ]                                   (0,1,0) to v
        //           [ u[2] v[2] w[2] ]                                   (0,0,1) to w

        A = R.transpose() * A;                                      //  A' = RotMat^T * A
        B = R.transpose() * B;                                      //  B' = RotMat^T * B
        C = R.transpose() * C;                                      //  C' = RotMat^T * C

        minX = min3(A(0), B(0), C(0));                              //  Find bounding box around translated, rotated triangle
        minY = min3(A(1), B(1), C(1));
        maxX = max3(A(0), B(0), C(0));
        maxY = max3(A(1), B(1), C(1));

        sampleLen = 0;
        p[1] = minY;
        while(p[1] <= maxY)
          {
            p[0] = minX;
            while(p[0] <= maxX)
              {
                q[0] = bary2d(C, B, p);
                q[1] = bary2d(B, A, p);
                q[2] = bary2d(A, C, p);

                if(q[0] >= 0.0 && q[1] >= 0.0 && q[2] >= 0.0)
                  {
                    if(++sampleLen == 1)                            //  Expand array
                      {
                        if((samples = (Vertex*)malloc(sizeof(Vertex))) == NULL)
                          {
                            cout << "ERROR: Unable to allocate sampled Vertex array." << endl;
                            return 0L;
                          }
                      }
                    else
                      {
                        if((samples = (Vertex*)realloc(samples, sampleLen * sizeof(Vertex))) == NULL)
                          {
                            cout << "ERROR: Unable to re-allocate sampled Vertex array." << endl;
                            return 0L;
                          }
                      }

                    samples[sampleLen - 1].x = (float)p[0];         //  Save (transformed) sample
                    samples[sampleLen - 1].y = (float)p[1];
                    samples[sampleLen - 1].z = 0.0;
                  }

                p[0] += sampleRate;
              }

            p[1] += sampleRate;
          }

        if(sampleLen > 0)                                           //  We have samples from this face to contribute to the running total
          {
            if(len == 0L)                                           //  Main total is so far empty
              {
                if(sampleLen >= SAMPLE_MAIN_MALLOC_MAX)             //  This contribution alone puts us over the limit
                  {
                    cacheSamples(points, len, cacheCount);          //  Write current accumulation to cache
                    len = 0L;                                       //  Reset the length of the main accumulator
                    free((*points));                                //  Release the accumulator array
                    if(((*points) = (Vertex*)malloc(sampleLen * sizeof(Vertex))) == NULL)
                      {
                        cout << "ERROR: Unable to allocate post-cache main sample array." << endl;
                        free(samples);
                        return 0L;
                      }
                  }
                else                                                //  This contribution (to zero) is within limit
                  {
                    if(((*points) = (Vertex*)malloc(sampleLen * sizeof(Vertex))) == NULL)
                      {
                        cout << "ERROR: Unable to allocate main sampled vertex array." << endl;
                        return 0L;
                      }
                  }
              }
            else                                                    //  Main total has already accumulated some samples
              {
                                                                    //  This contribution puts us over the limit
                if(len + (unsigned long)sampleLen >= SAMPLE_MAIN_MALLOC_MAX)
                  {
                    cacheSamples(points, len, cacheCount);          //  Write current accumulation to cache
                    len = 0L;                                       //  Reset the length of the main accumulator
                    free((*points));                                //  Release the accumulator array
                    if(((*points) = (Vertex*)malloc(sampleLen * sizeof(Vertex))) == NULL)
                      {
                        cout << "ERROR: Unable to allocate post-cache main sample array." << endl;
                        free(samples);
                        return 0L;
                      }
                  }
                else                                                //  This contribution is within limit
                  {
                    if(((*points) = (Vertex*)realloc((*points), (len + (unsigned long)sampleLen) * sizeof(Vertex))) == NULL)
                      {
                        cout << "ERROR: Unable to re-allocate main sampled vertex array." << endl;
                        return 0L;
                      }
                  }
              }

            for(j = 0; j < sampleLen; j++)                          //  Transform all in 'samples' back and save them to (*points)
              {
                A(0) = (double)samples[j].x;                        //  Save the sampled point in A
                A(1) = (double)samples[j].y;
                A(2) = (double)samples[j].z;
                //  RotMat^T = [ u[0] u[1] u[2] ]                       Takes u to (1,0,0)
                //             [ v[0] v[1] v[2] ]                             v to (0,1,0)
                //             [ w[0] w[1] w[2] ]                             w to (0,0,1)
                A = R * A;                                          //  A = RotMat * A'
                                                                    //  Restore the offset, vertex[0]
                (*points)[len + (unsigned long)j].x = A(0) + V[ F[i].v[0] ].x;
                (*points)[len + (unsigned long)j].y = A(1) + V[ F[i].v[0] ].y;
                (*points)[len + (unsigned long)j].z = A(2) + V[ F[i].v[0] ].z;
                (*points)[len + (unsigned long)j].r = V[ F[i].v[0] ].r;
                (*points)[len + (unsigned long)j].g = V[ F[i].v[0] ].g;
                (*points)[len + (unsigned long)j].b = V[ F[i].v[0] ].b;
                (*points)[len + (unsigned long)j].cam = false;
              }

            len += (unsigned long)sampleLen;                        //  Update length of target array
            free(samples);
          }
      }                                                             //  end for each Face

    if(len == 0L)
      cout << "Sampling completed without error, but no points were sampled." << endl;

    return len;
  }

/* At a fine enough sample rate, we will allocate all that there is to allocate.
   Call this function to write all samples to a cache file. */
void cacheSamples(Vertex** points, unsigned long len, unsigned int* cacheCtr)
  {
    ofstream fh;                                                    //  Output file handle
    char buffer[64];
    unsigned long i;
    float x, y, z;
    unsigned char r, g, b;

    sprintf(buffer, "psample.%d.cache", (*cacheCtr));

    #ifdef __PSAMPLE_DEBUG
    cout << "cacheSamples(" << buffer << ")" << endl;
    #endif

    fh.open(buffer);                                                //  Open the file
    fh.write(reinterpret_cast<char*>(&len), sizeof(long));          //  Write number of vertices in this cache file

    for(i = 0L; i < len; i++)                                       //  Write every Vertex
      {
        x = (*points)[i].x;
        y = (*points)[i].y;
        z = (*points)[i].z;

        r = (*points)[i].r;
        g = (*points)[i].g;
        b = (*points)[i].b;

        fh.write(reinterpret_cast<char*>(&x), sizeof(float));
        fh.write(reinterpret_cast<char*>(&y), sizeof(float));
        fh.write(reinterpret_cast<char*>(&z), sizeof(float));
        fh.write(reinterpret_cast<char*>(&r), sizeof(char));
        fh.write(reinterpret_cast<char*>(&g), sizeof(char));
        fh.write(reinterpret_cast<char*>(&b), sizeof(char));
      }

    fh.close();                                                     //  Close the file. Done!

    (*cacheCtr)++;                                                  //  Increase cache count

    return;
  }

/* Read the binary data of the given PLY file and load the given arrays.
   Return true on success.
   Return false on failure. */
bool loadMesh(const char* filename, unsigned int* info, Vertex** V, Face** F)
  {
    ifstream fh;                                                    //  File handle

    float x, y, z;                                                  //  Elements of Vertex
    unsigned char r, g, b;
    unsigned char isCam;
    unsigned char num;                                              //  Elements of Face
    unsigned int v0, v1, v2;
    unsigned char phase;                                            //  Into which array are we reading
    unsigned int i;
    bool breakout;

    #ifdef __PSAMPLE_DEBUG
    cout << "load(" << filename;
    cout << ", [" << +info[0] << ", " << +info[1] << ", " << +info[2] << ", " << +info[3] << "])" << endl;
    #endif

    fh.open(filename, ios::binary);                                 //  Open the given file (for binary reading this time)
    fh.seekg(info[3], ios::beg);                                    //  Jump past the header to the PLY data
                                                                    //  Allocate arrays for vertices, faces, tetrahedra
    if(((*V) = (Vertex*)malloc(info[0] * sizeof(Vertex))) == NULL)
      {
        cout << "ERROR: Unable to allocate vertex buffer." << endl;
        return false;
      }
    if(((*F) = (Face*)malloc(info[1] * sizeof(Face))) == NULL)
      {
        cout << "ERROR: Unable to allocate face buffer." << endl;
        free((*V));
        return false;
      }

    i = 0;
    phase = 0;
    breakout = false;
    while(!fh.eof() && !breakout)
      {
        switch(phase)
          {
            case 0:                                                 //  Read elements from file
                     fh.read(reinterpret_cast<char*>(&x), sizeof(float));
                     fh.read(reinterpret_cast<char*>(&y), sizeof(float));
                     fh.read(reinterpret_cast<char*>(&z), sizeof(float));
                     fh.read(reinterpret_cast<char*>(&r), sizeof(char));
                     fh.read(reinterpret_cast<char*>(&g), sizeof(char));
                     fh.read(reinterpret_cast<char*>(&b), sizeof(char));
                     fh.read(reinterpret_cast<char*>(&isCam), sizeof(char));
                     (*V)[i].x = x;                                 //  Write elements to struct in array
                     (*V)[i].y = y;
                     (*V)[i].z = z;
                     (*V)[i].r = r;
                     (*V)[i].g = g;
                     (*V)[i].b = b;
                     (*V)[i].cam = (isCam != 0) ? true : false;
                     #ifdef __PSAMPLE_LOAD_DEBUG
                     cout << "Vertex[" << +i << "]:" << endl;
                     cout << "    x: " << (*V)[i].x << " <-- " << x << endl;
                     cout << "    y: " << (*V)[i].y << " <-- " << y << endl;
                     cout << "    z: " << (*V)[i].z << " <-- " << z << endl;
                     cout << "    r: " << +(*V)[i].r << " <-- " << +r << endl;
                     cout << "    g: " << +(*V)[i].g << " <-- " << +g << endl;
                     cout << "    b: " << +(*V)[i].b << " <-- " << +b << endl;
                     if((*V)[i].cam)
                       cout << "  cam: 1";
                     else
                       cout << "  cam: 0";
                     cout << " <-- " << +isCam << endl;
                     #endif
                     if(++i == info[0])                             //  If we've read sufficient,
                       {                                            //  then change phases and reset the index
                         phase++;
                         i = 0;
                       }
                     break;
            case 1:                                                 //  Read elements from file
                                                                    //  Read 'num' just to advance the file pointer; we don't care.
                     fh.read(reinterpret_cast<char*>(&num), sizeof(char));
                     fh.read(reinterpret_cast<char*>(&v0), sizeof(int));
                     fh.read(reinterpret_cast<char*>(&v1), sizeof(int));
                     fh.read(reinterpret_cast<char*>(&v2), sizeof(int));
                     (*F)[i].v[0] = v0;
                     (*F)[i].v[1] = v1;
                     (*F)[i].v[2] = v2;
                     computeFaceBarycenter((*V), F, i);             //  Compute and store the barycenter for this Face
                     computeFaceNormal((*V), F, i);                 //  Compute and store the unit normal for this face
                     computeFaceArea((*V), F, i);                   //  Compute and store the area of this Face

                     #ifdef __PSAMPLE_LOAD_DEBUG
                     cout << "Face:" << endl;
                     cout << "  v0:     " << +(*F)[i].v[0] << " <-- " << +v0 << endl;
                     cout << "  v1:     " << +(*F)[i].v[1] << " <-- " << +v1 << endl;
                     cout << "  v2:     " << +(*F)[i].v[2] << " <-- " << +v2 << endl;
                     cout << "  normal: " << (*F)[i].normal[0] << ", " << (*F)[i].normal[1] << ", " << (*F)[i].normal[2] << endl;
                     cout << "  area:   " << (*F)[i].area << endl;
                     #endif
                     if(++i == info[1])                             //  If we've read sufficient,
                       breakout = true;                             //  then we're done!
                     break;
          }
      }
    fh.close();                                                     //  Close the file

    return true;
  }

/* Read the given PLY file and find out all the information we need to know about its binary data.
   Write the following details into the given array, 'info'.
     [0] Number of vertices                                         e.g.  element vertex 31480
     [1] Number of faces                                            e.g.  element face 408330
     [2] Number of tetrahedra                                       e.g.  element tetrahedra 204129
     [3] Offset of the first byte of PLY binary data
   Return true on success.
   Return false on failure. */
bool getMeshInfo(const char* filename, unsigned int** info)
  {
    ifstream fh;                                                    //  File handle
    unsigned int contentsLen;                                       //  Length of given file in bytes
    char* buffer;                                                   //  Contains file contents
    char line[256];                                                 //  Line buffer
    unsigned int i, j;
    unsigned int linectr;
    unsigned char spctr;                                            //  Count spaces
    char* word0;                                                    //  Substring: first word found in line
    char* word1;                                                    //  Substring: second word found in line
    char* word2;                                                    //  Substring: third word found in line
    unsigned char len0, len1, len2;                                 //  Lengths of these words

    #ifdef __PSAMPLE_DEBUG
    cout << "getMeshInfo(" << filename << ")" << endl;
    #endif

    fh.open(filename);                                              //  Open the given file
    fh.seekg(0, ios::end);                                          //  Jump to the end
    contentsLen = fh.tellg();                                       //  Find out how many bytes we need to allocate
    fh.seekg(0, ios::beg);                                          //  Jump back to the head of the file
    if((buffer = (char*)malloc(contentsLen * sizeof(char))) == NULL)//  Allocate the bytes
      return 0;
    fh.read(buffer, contentsLen);                                   //  Read in entire file
    fh.close();                                                     //  Close the file

    i = 0;                                                          //  We're only interested in the ASCII header for now,
    linectr = 0;                                                    //  so read line by line until "end_header" is found.
    while(i < contentsLen)                                          //  (This is an indulgent bound, but all we have so far.)
      {
        if(buffer[i] == '\n')                                       //  End of a line: time to examine the line buffer
          {
            line[linectr] = '\0';                                   //  Null-terminate the line buffer

            spctr = 0;                                              //  How many spaces on this header line?
            for(j = 0; j < linectr; j++)                            //  Count spaces
              {
                if(line[j] == ' ')
                  spctr++;
              }
            if(spctr == 2)                                          //  Two spaces found means three words on this line: a candidate!
              {
                j = 0;
                spctr = 0;                                          //  Reset space counter
                len0 = 0;                                           //  Set all initial lengths to zero
                len1 = 0;
                len2 = 0;
                while(j <= linectr)
                  {
                    switch(spctr)                                   //  Every iteration adds a character
                      {                                             //  to one of our three words.
                        case 0:  if(++len0 == 1)                    //  Expand the array
                                   {
                                     if((word0 = (char*)malloc(sizeof(char))) == NULL)
                                       return false;
                                   }
                                 else
                                   {
                                     if((word0 = (char*)realloc(word0, len0 * sizeof(char))) == NULL)
                                       return false;
                                   }
                                 //  Don't add the character just yet: check to see whether it's a space or a NULL!
                                 break;
                        case 1:  if(++len1 == 1)                    //  Expand the array
                                   {
                                     if((word1 = (char*)malloc(sizeof(char))) == NULL)
                                       return false;
                                   }
                                 else
                                   {
                                     if((word1 = (char*)realloc(word1, len1 * sizeof(char))) == NULL)
                                       return false;
                                   }
                                 //  Don't add the character just yet: check to see whether it's a space or a NULL!
                                 break;
                        default: if(++len2 == 1)                    //  Expand the array
                                   {
                                     if((word2 = (char*)malloc(sizeof(char))) == NULL)
                                       return false;
                                   }
                                 else
                                   {
                                     if((word2 = (char*)realloc(word2, len2 * sizeof(char))) == NULL)
                                       return false;
                                   }
                                 //  Don't add the character just yet: check to see whether it's a space or a NULL!
                      }

                    if(line[j] == ' ')                              //  Found a space: one word ends; point to the next.
                      {
                        switch(spctr)
                          {
                            case 0:  word0[len0 - 1] = '\0';  break;//  Cap word-0
                            case 1:  word1[len1 - 1] = '\0';  break;//  Cap word-1
                          }
                        spctr++;                                    //  Increment space counter
                      }
                    else if(line[j] == '\0')                        //  Found a null
                      word2[len2 - 1] = '\0';                       //  Cap word-2
                    else                                            //  Found a character
                      {                                             //  The current array has already been expanded,
                        switch(spctr)                               //  So just add the current character to the current array!
                          {
                            case 0:  word0[len0 - 1] = line[j];  break;
                            case 1:  word1[len1 - 1] = line[j];  break;
                            default: word2[len2 - 1] = line[j];  break;
                          }
                      }
                    j++;
                  }
                                                                    //  Now, is this actually what we want?
                if(strcmp(word0, "element") == 0 && strcmp(word1, "vertex") == 0)
                  (*info)[0] = (unsigned int)atoi(word2);
                if(strcmp(word0, "element") == 0 && strcmp(word1, "face") == 0)
                  (*info)[1] = (unsigned int)atoi(word2);
                if(strcmp(word0, "element") == 0 && strcmp(word1, "tetrahedra") == 0)
                  (*info)[2] = (unsigned int)atoi(word2);
              }

            linectr = 0;                                            //  Reset the line buffer's index, prepare for the next line

            if(strcmp(line, "end_header") == 0)                     //  We're done!
              {
                if(contentsLen > 0)
                  free(buffer);                                     //  Clean up
                (*info)[3] = ++i;                                   //  Point to the first byte of PLY data
                return true;                                        //  Let's get out of here!
              }
          }
        else
          {
            line[linectr] = buffer[i];                              //  Add this character to the line buffer
            linectr++;                                              //  Point to the next character
          }
        i++;
      }

    if(contentsLen > 0)
      free(buffer);                                                 //  Clean up

    return false;
  }

/* Compute and store the barycenter for this face. */
void computeFaceBarycenter(Vertex* V, Face** F, unsigned int f)
  {
    (*F)[f].barycenter[0] = (V[(*F)[f].v[0]].x + V[(*F)[f].v[1]].x + V[(*F)[f].v[2]].x) / 3.0;
    (*F)[f].barycenter[1] = (V[(*F)[f].v[0]].y + V[(*F)[f].v[1]].y + V[(*F)[f].v[2]].y) / 3.0;
    (*F)[f].barycenter[2] = (V[(*F)[f].v[0]].z + V[(*F)[f].v[1]].z + V[(*F)[f].v[2]].z) / 3.0;
    return;
  }

/* Area of a 3D triangle is 1/2 times the magnitude of the cross product of two of its vectors. */
void computeFaceArea(Vertex* V, Face** F, unsigned int f)
  {
    float* u;
    float* v;
    float* cross;

    if((u = (float*)malloc(3 * sizeof(float))) == NULL)
      {
        cout << "ERROR: Unable to allocate cross-product input array u." << endl;
        exit(1);
      }
    if((v = (float*)malloc(3 * sizeof(float))) == NULL)
      {
        cout << "ERROR: Unable to allocate cross-product input array v." << endl;
        exit(1);
      }
    if((cross = (float*)malloc(3 * sizeof(float))) == NULL)
      {
        cout << "ERROR: Unable to allocate cross-product output array." << endl;
        exit(1);
      }

    u[0] = V[(*F)[f].v[0]].x - V[(*F)[f].v[1]].x;
    u[1] = V[(*F)[f].v[0]].y - V[(*F)[f].v[1]].y;
    u[2] = V[(*F)[f].v[0]].z - V[(*F)[f].v[1]].z;

    v[0] = V[(*F)[f].v[0]].x - V[(*F)[f].v[2]].x;
    v[1] = V[(*F)[f].v[0]].y - V[(*F)[f].v[2]].y;
    v[2] = V[(*F)[f].v[0]].z - V[(*F)[f].v[2]].z;

    crossProduct(u, v, &cross);

    (*F)[f].area = 0.5 * sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);

    free(u);
    free(v);
    free(cross);

    return;
  }

/* Compute and store the unit normal for the Face at index 'f' in array 'F'. */
void computeFaceNormal(Vertex* V, Face** F, unsigned int f)
  {
    double Nx, Ny, Nz;
    double mag;

    //  (B - A) cross (C - A)

    //  [a]       [d]     [bf - ce]
    //  [b] cross [e]  =  [cd - af]
    //  [c]       [f]     [ae - bd]

    //             B                   A
    //  a = V[(*F)[f].v[1]].x - V[(*F)[f].v[0]].x
    //  b = V[(*F)[f].v[1]].y - V[(*F)[f].v[0]].y
    //  c = V[(*F)[f].v[1]].z - V[(*F)[f].v[0]].z

    //             C                   A
    //  d = V[(*F)[f].v[2]].x - V[(*F)[f].v[0]].x
    //  e = V[(*F)[f].v[2]].y - V[(*F)[f].v[0]].y
    //  f = V[(*F)[f].v[2]].z - V[(*F)[f].v[0]].z

    Nx = (V[(*F)[f].v[1]].y - V[(*F)[f].v[0]].y) * (V[(*F)[f].v[2]].z - V[(*F)[f].v[0]].z) - (V[(*F)[f].v[1]].z - V[(*F)[f].v[0]].z) * (V[(*F)[f].v[2]].y - V[(*F)[f].v[0]].y);
    Ny = (V[(*F)[f].v[1]].z - V[(*F)[f].v[0]].z) * (V[(*F)[f].v[2]].x - V[(*F)[f].v[0]].x) - (V[(*F)[f].v[1]].x - V[(*F)[f].v[0]].x) * (V[(*F)[f].v[2]].z - V[(*F)[f].v[0]].z);
    Nz = (V[(*F)[f].v[1]].x - V[(*F)[f].v[0]].x) * (V[(*F)[f].v[2]].y - V[(*F)[f].v[0]].y) - (V[(*F)[f].v[1]].y - V[(*F)[f].v[0]].y) * (V[(*F)[f].v[2]].x - V[(*F)[f].v[0]].x);
    mag = sqrt(Nx * Nx + Ny * Ny + Nz * Nz);

    (*F)[f].normal[0] = Nx / mag;
    (*F)[f].normal[1] = Ny / mag;
    (*F)[f].normal[2] = Nz / mag;

    return;
  }

/* Compute the cross-product, stick the result in the given array 'output' and return the magnitude. */
float crossProduct(float* u, float* v, float** output)
  {
    (*output)[0] = u[1] * v[2] - u[2] * v[1];
    (*output)[1] = u[2] * v[0] - u[0] * v[2];
    (*output)[2] = u[0] * v[1] - u[1] * v[0];
    return sqrt((*output)[0] * (*output)[0] + (*output)[1] * (*output)[1] + (*output)[2] * (*output)[2]);
  }

/* Determine 2D barycentric coordinates for point C */
double bary2d(Vector3f A, Vector3f B, double* C)
  {
    return (B(0) - A(0)) * (C[1] - A(1)) - (B(1) - A(1)) * (C[0] - A(0));
  }

/* Return the least of three floats */
double min3(double a, double b, double c)
  {
    double m = a;
    if(b < m)
      m = b;
    if(c < m)
      m = c;
    return m;
  }

/* Return the greatest of three floats */
double max3(double a, double b, double c)
  {
    double m = a;
    if(b > m)
      m = b;
    if(c > m)
      m = c;
    return m;
  }

void usage(void)
  {
    cout << "Usage:  ./psample ply-filename <sample-rate>" << endl;
    cout << " e.g.:  ./psample terrains.ply" << endl;
    cout << " e.g.:  ./psample terrains.ply 0.001" << endl;
    cout << "Writes a new file \"psample.ply\" containing only vertices sampled from the faces in the source." << endl;
    return;
  }
