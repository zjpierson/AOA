////////////////////////////////////////////////////////////////////
//                                                                //
//               Quicksort With Median of 3 Pivot                 //
//                                                                //
//   Uses stack instead of recursive calls                        //
//   Makes use of tail recursion to simplify stack use            //
//   Processes small half first to minimize stack space used      //
//                                                                //
////////////////////////////////////////////////////////////////////

//#define DEBUG
#include <iostream>
#include <fstream>
#include <math.h>
#include <stack>
#include <cstdlib>
#include <values.h>  // for MAXLONG
#include <sys/resource.h>  // for setrlimit (stack size)
using namespace std;

const int  MAX_N =  1000000;
unsigned long comparisons = 0;

////////////////////////////////////////////////////////////////////

float clock_seconds()
  {
  return clock() / (float) CLOCKS_PER_SEC;
  }

////////////////////////////////////////////////////////////////////

//Generates a list that's already in order
void gen_rand_list(long a[], int n)
  {
  int  i;
  long x;

  x = 0;
  for (i = 1; i < n; i++)
    {
      x += rand() % 100;
      a[i] = x;
    }
  }

/////////////////////////////////////////////////////////////////////

void gen_worst_case_md3(long a[], int n)
{
    long long int i;

    for(i = 1; i < n-1; i++)
        a[i] = i+1;
    a[0] = 0;
    a[n-1] = 1;


#ifdef DEBUG
    for(i = 0; i < n; i++)
        cout << a[i] << "  ";
    cout << endl;
#endif
}

////////////////////////////////////////////////////////////////////

void check_order(long a[], int n)
  {
  int i;
  bool ordered;

  ordered = true;
  for(i = 0; i < n - 1; i++)
    {
    if (a[i] > a[i + 1])
      ordered = false;
    }
  if (ordered)
    cout << "List is in order" << endl;
  else
    cout << "List is NOT in order" << endl;
  }



////////////////////////////////////////////////////////////////////

void swap(long &x, long &y)
  {
  long oldx;
  long oldy;

  oldx = x;
  oldy = y;
  x = oldy;
  y = oldx;
  }

////////////////////////////////////////////////////////////////////

void print_list(long a[], int n)
  {
  int i;

  for(i = 0; i < n; i++)
    {
    cout << " " << a[i];
    if ((i + 1) % 10 == 0)
      cout << endl;
    }
  cout << endl;
  }

////////////////////////////////////////////////////////////////////
//   swap median of first, middle, last to first so               //
//   quicksort function does not need to be changed               //
////////////////////////////////////////////////////////////////////

void median_of_3_to_left(long a[], int left, int right)
  {
  long a_left;
  long a_right;
  long a_mid;
  int  mid;

  mid = (left + right) / 2;
  a_left = a[left];
  a_right = a[right];
  a_mid = a[mid];

  // check for left value already median
  comparisons++;
  if (a_left >= a_mid)
  {
    comparisons++;
    if (a_left <= a_right)
      return;
  }

  comparisons++;
  if (a_left <= a_mid)
  {
    comparisons++;
    if (a_left >= a_right)
      return;
  }

  // check for middle value is median
  comparisons++;
  if (a_mid >= a_left)
  {
    comparisons++;
    if (a_mid <= a_right)
    {
      swap(a[mid], a[left]);
      return;
    }
  }

  comparisons++;
  if (a_mid <= a_left)
  {
    comparisons++;
    if (a_mid >= a_right)
    {
      swap(a[mid], a[left]);
      return;
    }
  }

  // others aren't, so right value must be median
  swap(a[right], a[left]);
  }

////////////////////////////////////////////////////////////////////
//                          quicksort                             //
//    uses original partition algorithm                           //
//    uses median of three as pivot                               //
////////////////////////////////////////////////////////////////////

void quick_sort_median_3(long a[], int left, int right)
  {
  int  i;
  int  j;
  long pivot;

  comparisons++;
  if (left < right)
    {
    i = left;
    j = right + 1;
    median_of_3_to_left(a, left, right);
    pivot = a[left];

    do
    {
      comparisons++;
      do
      {
        comparisons++;
        i++;
      }
      while(a[i] < pivot);

      do
      {
        comparisons++;
        j--;
      }
      while(a[j] > pivot);

      comparisons++;
      if (i < j)
        swap(a[i], a[j]);

    } // end do
    while (i <= j);
    

    swap(a[left], a[j]);
    quick_sort_median_3(a, left, j - 1);
    quick_sort_median_3(a, j + 1, right);
    }  // end if

  } // quick sort

////////////////////////////////////////////////////////////////////

struct stack_frame
  {
  int  left;
  int  right;
  };

////////////////////////////////////////////////////////////////////

void quick_sort_left(long a[], int left, int right)
  {
  int  i;
  int  j;
  long pivot;
  stack <stack_frame> s;
  stack_frame current;
  stack_frame left_side;
  stack_frame right_side;

  current . left = left;
  current . right = right;
  s . push(current);

  comparisons++;
  while (!s . empty())
    {
    comparisons++;
    current = s . top();
    s . pop();
    left = current . left;
    right = current . right;

    comparisons++;
    if (left < right)
      {
      i = left;
      j = right + 1;
      pivot = a[left];

      do
      {
        comparisons++;
        do
        {
          comparisons++;
          i++;
        }
        while(a[i] < pivot);

        do
        {
          comparisons++;
          j--;
        }
        while(a[j] > pivot);

        comparisons++;
        if (i < j)
          swap(a[i], a[j]);

      } // end do
      while (i <= j);

      swap(a[left], a[j]);

      // stack bigger half and process smaller half next
      // quick_sort_left(a, left, j - 1);
      left_side . left = left;
      left_side . right = j - 1;
      // quick_sort_left(a, j + 1, right);
      right_side . left = j + 1;
      right_side . right = right;

      comparisons++;
      if ((j - 1) - left + 1 > right - (j + 1) + 1)
        {
        // left side bigger
        s . push(left_side);
        s . push(right_side);
        }
      else
        {
        // right side bigger
        s . push(right_side);
        s . push(left_side);
        }
      }  // end if
    } // while stack not empty

  } // quick sort

////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
  {
  long      a[MAX_N + 1]; // quicksort needs an extra spot at end
  long      b[MAX_N + 1]; // copy for sorting again
  int       i;
  int       n;
  float     t1;
  float     t2;
  float     t3;
  float     t4;
  rlimit    rlim;
  ofstream  times;
  ofstream  bad;
/*********** Set stack size as large as possible  ***
****************/
  rlim . rlim_cur = RLIM_INFINITY;
  rlim . rlim_max = RLIM_INFINITY;
/********* Set small stack for testing stack overflow **
  rlim . rlim_cur = 24000000;
  rlim . rlim_max = 24000000;
*/
  setrlimit(RLIMIT_STACK,&rlim);
  
  //check cmd line args
  if(argc == 1)
    n = 200000;
  else if( argc == 2)
    n = atoi(argv[1]);
  else
  {
    cout << "USAGE:\tquickSort  [length of list]i" << endl;
    return 1;
  }


  cout << "Running quicksorts with " << n << " numbers" << endl;
  if (rlim . rlim_max == RLIM_INFINITY)
    cout << "Using maximum stack size." << endl;
  else
    cout << "Using stack size " << rlim . rlim_max << "." << endl;
  
  //open output files
  times.open("quicksort.times");
  bad.open("quicksort.bad");

//    generate list in order
//    gen_rand_list(a, n);
   gen_worst_case_md3(a, n);

   //print worst case list to file "quicksort.bad"
   for(i = 0; i < n; i++)
   {
       bad << a[i] << " "; 
   }

  //copy list a into list b
  for (i = 0; i < n; i++)
    b[i] = a[i];

  //set end of list to MAXLONG
  a[n] = MAXLONG;
  b[n] = MAXLONG;

  cout << "Before quicksort_median_3" << endl;
  t1 = clock_seconds();
  quick_sort_median_3(a, 0, n - 1);
  t2 = clock_seconds();
  cout << "After quicksort_median_3" << endl;
  check_order(a, n);
  cerr << "median 3 took " << t2 - t1 << " seconds." << endl;

  cout << "Comparisons: " << comparisons << endl;
  times << "Median of 3 quicksort Comparisons: " << comparisons << endl;
  times << "Median of 3 quicksort time: " << t2-t1 << " seconds." << endl;
  comparisons = 0;

  cout << "Before quicksort_left" << endl;
  a[n] = MAXINT;
  t3 = clock_seconds();
  quick_sort_left(b, 0, n - 1);
  t4 = clock_seconds();
  cout << "After quicksort_left" << endl;
  check_order(b, n);
  cerr << "left took " << t4 - t3 << " seconds." << endl;

  cout << "Comparisons: " << comparisons << endl;
  times << "Quicksort_left Comparisons: " << comparisons << endl;
  times << "Quicksort_left times: " << t4 - t3 << " seconds." << endl;
  cout << "Done with quicksort" << endl;

  bad.close();
  times.close();

  return 0;
  }

////////////////////////////////////////////////////////////////////
//                End median of 3 quicksort                       //
////////////////////////////////////////////////////////////////////

