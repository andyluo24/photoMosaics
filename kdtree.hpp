/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    /**
     * @todo Implement this function!
     */
     if (first[curDim] < second[curDim]) {
       return true;
     } else if (first[curDim] > second[curDim]) {
       return false;
     } else {
       return first < second;
     }
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    /**
     * @todo Implement this function!
     */
     double d_pt1 = 0;
     for (unsigned i = 0; i < Dim; i++) {
       d_pt1 += pow(potential[i]-target[i], 2);
     }
     double d_pt = sqrt(d_pt1);
     double d_ct1 = 0;
     for (unsigned j = 0; j < Dim; j++) {
       d_ct1 += pow(currentBest[j]-target[j], 2);
     }
     double d_ct = sqrt(d_ct1);
     if (d_pt < d_ct) {
       return true;
     } else if (d_pt > d_ct) {
       return false;
     } else {
       return potential < currentBest;
     }
}

template <int Dim>
void KDTree<Dim>::construct(int start, int end, int dimension, KDTreeNode *& subRoot)
{
  int median_index = (start + end) / 2;
  if (start > end) {
    return;
  }
  quickSelect(start, end, median_index, dimension);
  subRoot = new KDTreeNode(points[median_index]);
  construct(start, median_index - 1, (dimension + 1) % Dim, subRoot->left);
  construct(median_index + 1, end, (dimension + 1) % Dim, subRoot->right);
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
    /**
     * @todo Implement this function!
     */
   size = newPoints.size();
   points = newPoints;
   root = NULL;
   if(size == 0) {
     return;
   }
   construct(0, size - 1, 0, root);
}

template <int Dim>
typename KDTree<Dim>::KDTreeNode *KDTree<Dim>::copyHelper(const KDTreeNode *& node)
{
  if (node != NULL) {
     KDTreeNode * newNode = new KDTreeNode(node->point);
     newNode->left = copyHelper(node->left);
     newNode->right = copyHelper(node->right);
     return newNode;
   }
   return NULL;
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  /**
   * @todo Implement this function!
   */
  points = other.points;
  size = other.size;
  copyHelper(other.root);
}

template <int Dim>
int KDTree<Dim>::partition(int start, int end, int dimension)
{
  Point<Dim> pivot = points[end];
  int i = start;
  for (int j = start; j < end; j++) {
    if (smallerDimVal(points[j], pivot, dimension)) {
      swap(points[i], points[j]);
      i++;
    }
  }
  swap(points[i], points[end]);
  return i;
}

template <int Dim>
void KDTree<Dim>::quickSelect(int start, int end, int k, int dimension)
{
  if (start >= end) {
    return;
  }
  int pivot_index = partition(start, end, dimension);
  if (k < pivot_index) {
    quickSelect(start, pivot_index - 1, k, dimension);
  }
  if (k > pivot_index) {
    quickSelect(pivot_index + 1, end, k, dimension);
  }
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  /**
   * @todo Implement this function!
   */
  if (this != &rhs) {
    clearHelper();
    copyHelper(rhs);
  }
  return *this;
}

template <int Dim>
void KDTree<Dim>::clearHelper(KDTreeNode *& root)
{
   if (root != NULL) {
     if (root->left != NULL) {
       clearHelper(root->left);
     }
     if (root->right != NULL) {
       clearHelper(root->right);
     }
     delete root;
     root = NULL;
   }
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  /**
   * @todo Implement this function!
   */
   clearHelper(root);
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */
    return findHelper(root, 0, query);
}
template <int Dim>
Point<Dim> KDTree<Dim>::findHelper(KDTreeNode * r, int d, const Point<Dim>& target) const {
  Point<Dim> checkup = r->point;
  if (r->left == NULL && r->right == NULL) {
    return checkup;
  }

  Point<Dim> traveler = r->point;
	bool check;

	if (!smallerDimVal(target, r->point, d)) {
    traveler = (r->right == nullptr) ? findHelper(r->left, (d + 1) % Dim, target) : findHelper(r->right, (d + 1) % Dim, target);
	} else {
    traveler = (r->left == nullptr) ? findHelper(r->right, (d + 1) % Dim, target) : findHelper(r->left, (d + 1) % Dim, target);
	}

  check = (!smallerDimVal(target, r->point, d)) ? false : true;

	if (shouldReplace(target, traveler, r->point) == true) {
    traveler = r->point;
	}

	double banjing = 0;
  unsigned j = 0;
  while (j < Dim) {
    double banjing1 = pow((target[j] - traveler[j]), 2);
    banjing += banjing1;
    j++;
  }

	double cuttingedge = r->point[d] - target[d];
	cuttingedge = pow(cuttingedge, 2);

	if (cuttingedge < banjing || cuttingedge == banjing) {
    KDTreeNode* c = nullptr;
    if (check == true) {
      c = r->right;
    } else {
      c = r->left;
    }

		if (c != NULL) {
      Point<Dim> *backup = nullptr;
			Point<Dim> ano = findHelper(c, (d + 1) % Dim, target);
      backup = &ano;
		  if (shouldReplace(target, traveler, *backup) == true) {
        traveler = *backup;
      }
		}
	}
	return traveler;
}
