/* Eric Joyce, Stevens Institute of Technology, 2019; ejoyce@stevens.edu

   19dec19:  Smoothes the given mesh vertices.
             Assumes vertices have attributes x, y, z, r, g, b, and isCam
*/

#include <fstream>
#include <iostream>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
#define __SMOOTH_DEBUG 1
*/
/*
#define __SMOOTH_LOAD_DEBUG 1
*/

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

typedef struct VertexLookupType
  {
    unsigned int* f;                                                //  Array of Faces (all Faces containing Vertex matching this struct's index
    unsigned int len;                                               //  Length of that array
    unsigned int ptr;                                               //  Which cell of the array we're pointing at
  } VertexLookup;

/*
g++ -Wall smooth.cpp -lm -o smooth
./smooth terrains.mincut.KB.ply 0.5
*/

using namespace std;

bool smooth(Vertex**, Face*, unsigned int*, double);
void writePLY(const char*, Vertex*, Face*, unsigned int*);
bool loadMesh(const char*, unsigned int*, Vertex**, Face**);
bool getMeshInfo(const char*, unsigned int**);
void computeFaceBarycenter(Vertex*, Face**, unsigned int);
void computeFaceArea(Vertex*, Face**, unsigned int);
void computeFaceNormal(Vertex*, Face**, unsigned int);
float crossProduct(float*, float*, float**);
void usage(void);

int main(int argc, char* argv[])
  {
    Vertex* V;                                                      //  Array of vertices
    Face* F;                                                        //  Array of faces

    double smoothRate = 0.5;                                        //  Scalar by which we tune the smoothing factor.
    unsigned int* fileInfo;                                         //  Array of data about the given PLY file
    unsigned char fileInfoLen;                                      //  Length of that array

    fileInfoLen = 4;                                                //  There are 4 things we need to know about the PLY:
    if(argc < 2)                                                    //  [0] Number of vertices
      {                                                             //  [1] Number of faces
        usage();                                                    //  [2] Number of tetrahedra (not needed here)
        return 0;                                                   //  [3] Offset of the first byte of PLY binary data
      }                                                             //  And make sure we have enough arguments!
    if(argc > 2)                                                    //  This argument is optional:
      smoothRate = (double)atof(argv[2]);                           //  If given a 'smoothRate', use that

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
        if(fileInfo[0] > 0)                                         //  Dump 'n' croak
          free(V);
        if(fileInfo[1] > 0)
          free(F);
        free(fileInfo);
        return 0;
      }

    if(!smooth(&V, F, fileInfo, smoothRate))                        //  Smoooooooooth
      {
        if(fileInfo[0] > 0)                                         //  Dump 'n' croak
          free(V);
        if(fileInfo[1] > 0)
          free(F);
        free(fileInfo);
        return 0;
      }

    writePLY("smooth.ply", V, F, fileInfo);                         //  Create the new file

    if(fileInfo[0] > 0)                                             //  These are no longer necessary
      free(V);
    if(fileInfo[1] > 0)
      free(F);
    free(fileInfo);

    return 0;
  }

/* Smooth all vertices in 'V' by factor 's'.
   'F' and 'info' are passed for reference only. */
bool smooth(Vertex** V, Face* F, unsigned int* info, double s)
  {
    Vertex* V2;                                                     //  New Vertex array
                                                                    //  Vertex Lookup Table:
                                                                    //  key = vertex #; val = all faces containing this vertex
    VertexLookup* lookup;                                           //  Length = info[0] = total number of Vertices
    unsigned int* neighbors;                                        //  Array of Vertex neighbors
    unsigned int nlen;                                              //  Length of that array
    unsigned int face;                                              //  Convenient Face index holder
    unsigned int i, j, k, n;
    float x, y, z;                                                  //  Average of neighbors
    float dx, dy, dz;                                               //  Vector from original point to average point

    #ifdef __SMOOTH_DEBUG
    cout << "smooth(" << s << ")" << endl;
    #endif

    if((lookup = (VertexLookup*)malloc(info[0] * sizeof(VertexLookup))) == NULL)
      {
        cout << "ERROR: Unable to allocate vertex lookup array." << endl;
        return false;
      }
    for(i = 0; i < info[0]; i++)                                    //  Initialize lookup table:
      {
        lookup[i].len = 0;                                          //  Before reading, all Vertices belong to zero Faces
        lookup[i].ptr = 0;                                          //  and their accumulators point to the beginning.
      }
    for(i = 0; i < info[1]; i++)                                    //  For every Face...
      {
        for(j = 0; j < 3; j++)                                      //  For every Face Vertex...
          lookup[ F[i].v[j] ].len++;                                //  increase the number of Faces in which Vertex appears
      }
    for(i = 0; i < info[0]; i++)                                    //  We have counts, now allocate arrays
      {
        if((lookup[i].f = (unsigned int*)malloc(lookup[i].len * sizeof(int))) == NULL)
          {
            cout << "ERROR: Unable to allocate array of member Faces for Vertex lookup[" << +i << "]." << endl;
            for(j = 0; j < i; j++)                                  //  Release arrays allocated previous to this failure
              {
                if(lookup[j].len > 0)
                  free(lookup[j].f);
              }
            free(lookup);                                           //  Release the main stash
            return false;                                           //  Slink home in shame
          }
      }
    for(i = 0; i < info[1]; i++)                                    //  FINALLY, build the actual lookup table:
      {                                                             //  For every Face...
        for(j = 0; j < 3; j++)                                      //  For every Vertex in that Face...
          {
            lookup[ F[i].v[j] ].f[ lookup[ F[i].v[j] ].ptr ] = i;   //  Add to the cell to which we currently point
            lookup[ F[i].v[j] ].ptr++;                              //  Move the pointer forward
          }
      }

    if((V2 = (Vertex*)malloc(info[0] * sizeof(Vertex))) == NULL)    //  Allocate a new array of SMOOTHED Vertices
      {
        cout << "ERROR: Unable to allocate array for smoothed Vertices." << endl;
        for(i = 0; i < info[0]; i++)                                //  Free lookup arrays
          {
            if(lookup[i].len > 0)
              free(lookup[i].f);
          }
        if(info[0] > 0)
          free(lookup);                                             //  Free lookup itself
        return false;                                               //  Slink home in shame
      }

    for(i = 0; i < info[0]; i++)                                    //  For each Vertex i in V...
      {
        nlen = 0;                                                   //  Find all Vertex i's neighbors:
        for(j = 0; j < lookup[i].len; j++)                          //  Go through each Face that contains Vertex i
          {
            face = lookup[i].f[j];                                  //  Save a convenient reference to that Face
            for(k = 0; k < 3; k++)                                  //  Go through all of that Face's Vertices
              {
                if(F[face].v[k] != i)                               //  Only consider Vertices other than Vertex i
                  {
                    n = 0;                                          //  Make sure this Vertex is not already listed as a neighbor
                    while(n < nlen && neighbors[n] != F[face].v[k])
                      n++;
                    if(n == nlen)                                   //  It is NOT already listed:
                      {                                             //  add it now
                        if(++nlen == 1)
                          {
                            if((neighbors = (unsigned int*)malloc(sizeof(int))) == NULL)
                              {
                                cout << "ERROR: Unable to allocate Vertex neighbors array." << endl;
                              }
                          }
                        else
                          {
                            if((neighbors = (unsigned int*)realloc(neighbors, nlen * sizeof(int))) == NULL)
                              {
                                cout << "ERROR: Unable to re-allocate Vertex neighbors array." << endl;
                              }
                          }
                        neighbors[nlen - 1] = F[face].v[k];
                      }
                  }
              }                                                     //  Done collecting neighbors from 'face'
          }                                                         //  Done collecting Vertex neighbors

        x = 0.0;  y = 0.0;  z = 0.0;                                //  Start a sum of neighbors
        for(j = 0; j < nlen; j++)                                   //  Add neighbors together
          {
            x += (*V)[ neighbors[j] ].x;
            y += (*V)[ neighbors[j] ].y;
            z += (*V)[ neighbors[j] ].z;
          }
        x /= (float)nlen;                                           //  Average neighbors together
        y /= (float)nlen;
        z /= (float)nlen;

        dx = x - (*V)[i].x;                                         //  Find vector
        dy = y - (*V)[i].y;
        dz = z - (*V)[i].z;

        V2[i].x = (*V)[i].x + dx * s;                               //  New point moves 's' along 'dx' toward 'x'
        V2[i].y = (*V)[i].y + dy * s;                               //  New point moves 's' along 'dy' toward 'y'
        V2[i].z = (*V)[i].z + dz * s;                               //  New point moves 's' along 'dz' toward 'z'

        if(nlen > 0)                                                //  Done with this Vertex;
          free(neighbors);                                          //  release the neighbors
      }

    for(i = 0; i < info[0]; i++)                                    //  Update values in array 'V' with smoothed values
      {
        (*V)[i].x = V2[i].x;
        (*V)[i].y = V2[i].y;
        (*V)[i].z = V2[i].z;
      }

    if(info[0] > 0)                                                 //  Clean up:
      {
        free(V2);                                                   //  Free smoothed tmp array

        for(i = 0; i < info[0]; i++)                                //  Free lookup arrays
          {
            if(lookup[i].len > 0)
              free(lookup[i].f);
          }
        free(lookup);                                               //  Free lookup itself
      }

    return true;
  }

/* Smoothing is complete, and we write the mesh to a new file.
   Notice that though our work in this program did not change the Faces, we still
   receive their array as an argument so that we can write them to file.
   Also note that here we DO include the isCam attribute. */
void writePLY(const char* filename, Vertex* V, Face* F, unsigned int* info)
  {
    ofstream fh;                                                    //  Output file handle
    unsigned int i, v0, v1, v2;
    unsigned char num = 3, boolHolder;

    #ifdef __SMOOTH_DEBUG
    cout << "writePLY()" << endl;
    #endif

    fh.open(filename);                                              //  Open the file
    fh << "ply" << endl;                                            //  Ply format header
    fh << "format binary_little_endian 1.0" << endl;                //  Declare endianness
    fh << "element vertex " << +info[0] << endl;                    //  Declare number of vertices
    fh << "property float x" << endl;                               //  Vertex properties are the same as in the original PLY
    fh << "property float y" << endl;
    fh << "property float z" << endl;
    fh << "property uchar red" << endl;
    fh << "property uchar green" << endl;
    fh << "property uchar blue" << endl;
    fh << "property uchar isCam" << endl;                           //  Include this because subsampling takes it out
    fh << "element face " << +info[1] << endl;
    fh << "property list uchar int vertex_index" << endl;           //  No labels, no photoconsistency here
    fh << "end_header" << endl;                                     //  Close the header

    for(i = 0; i < info[0]; i++)                                    //  Write every Vertex in 'V' array
      {
        fh.write(reinterpret_cast<char*>(&V[i].x), sizeof(float));
        fh.write(reinterpret_cast<char*>(&V[i].y), sizeof(float));
        fh.write(reinterpret_cast<char*>(&V[i].z), sizeof(float));
        fh.write(reinterpret_cast<char*>(&V[i].r), sizeof(char));
        fh.write(reinterpret_cast<char*>(&V[i].g), sizeof(char));
        fh.write(reinterpret_cast<char*>(&V[i].b), sizeof(char));
        boolHolder = (V[i].cam) ? 1 : 0;
        fh.write(reinterpret_cast<char*>(&boolHolder), sizeof(char));
      }

    for(i = 0; i < info[1]; i++)                                    //  Write every Face in 'F' array
      {
        v0 = F[i].v[0];
        v1 = F[i].v[1];
        v2 = F[i].v[2];
        fh.write(reinterpret_cast<char*>(&num), sizeof(char));
        fh.write(reinterpret_cast<char*>(&v0), sizeof(int));
        fh.write(reinterpret_cast<char*>(&v1), sizeof(int));
        fh.write(reinterpret_cast<char*>(&v2), sizeof(int));
      }

    fh.close();                                                     //  Close the file. Done!

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

    #ifdef __SMOOTH_DEBUG
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
                     #ifdef __SMOOTH_LOAD_DEBUG
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

                     #ifdef __SMOOTH_LOAD_DEBUG
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

    #ifdef __SMOOTH_DEBUG
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

void usage(void)
  {
    cout << "Usage:  ./smooth ply-filename <smoothing-rate>" << endl;
    cout << " e.g.:  ./smooth terrains.ply" << endl;
    cout << " e.g.:  ./smooth terrains.ply 0.01" << endl;
    cout << "Default smoothing rate is 0.5" << endl;
    cout << "Writes a new file \"smooth.ply\" with smoothed vertices." << endl;
    return;
  }
