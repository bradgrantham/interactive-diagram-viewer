#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace std;

const int max_node_rows = 100;

struct Node {
    int parented;
    string name;
    float x, y, weight;
    vector<Node*> sibs;
    Node *next;

    Node(const string& name_, Node* unplaced) :
        parented(0),
        name(name_),
        x(0.0f),
        y(0.0f),
        weight(0.0f),
        next(unplaced)
    {}
};

Node *nodesUnplaced;
int nodesUnplacedCount;
Node *nodeRows[max_node_rows];
int nodeRowHeights[max_node_rows];
int nodeRowCount;
int tallestRow, tallestRowHeight = -1;

Node *findNode(const string& name)
{
    Node *cur;

    cur = nodesUnplaced;
    while(cur != NULL) {
	if(cur->name == name)
	    return cur;
	cur = cur->next;
    }
    cur = new Node(name, nodesUnplaced);
    nodesUnplaced = cur;
    nodesUnplacedCount++;
    return cur;
}

void makeSib(Node *a, Node *b)
{
    for(auto it = a->sibs.begin(); it != a->sibs.end(); it++)
	if(*it == b)
	    return;
    a->sibs.push_back(b);
}

void readNodes()
{
    char word1[512], word2[512];
    while(scanf("%s %s", word1, word2) == 2) {
        Node *first = findNode(string(word1));
        Node *second = findNode(word2);
	
	makeSib(first, second);
	second->parented++;
    }
}

void orderNodesInRows()
{
    Node *cur, **prevPtr;
    Node *placedTmp;
    int row = 0, rowNodeCount;

    while(nodesUnplaced != NULL) {
	prevPtr = &nodesUnplaced;
	cur = nodesUnplaced;
	placedTmp = NULL;
	rowNodeCount = 0;
	while(cur != NULL) {
	    if(cur->parented == 0) {
		*prevPtr = cur->next;
		cur->next = nodeRows[row];
		nodeRows[row] = cur;
		cur = *prevPtr;
		rowNodeCount++;
		nodesUnplacedCount--;
	    } else {
		prevPtr = &cur->next;
		cur = cur->next;
	    }
	}

	if(rowNodeCount == 0 && nodesUnplacedCount > 0) {
	    fprintf(stderr, "Graph cycle encountered with the following"
		" nodes remaining:\n");
	    cur = nodesUnplaced;
	    while(cur != NULL) {
		fprintf(stderr, "    %s\n", cur->name.c_str());
		cur = cur->next;
	    }
	    fprintf(stderr, "These nodes will be discarded.\n");
	    nodesUnplaced = NULL;
	}

	cur = nodeRows[row];
	while(cur != NULL) {
            for(auto it = cur->sibs.begin(); it != cur->sibs.end(); it++)
		(*it)->parented--;

	    cur = cur->next;
	}

	nodeRowHeights[row] = rowNodeCount;

	if(rowNodeCount > tallestRowHeight) {
	    tallestRowHeight = rowNodeCount;
	    tallestRow = row;
	}
	row++;
    }
    nodeRowCount = row;
}

void placeNodes()
{
    Node *cur;
    float x, y, starty = 0.0f;

    if(true) {
        x = 0; 
        for(int i = 0; i < nodeRowCount; i++) {
            cur = nodeRows[i];
            y = (tallestRowHeight - nodeRowHeights[i]) / 2.0f;
            while(cur != NULL) {
                cur->x = x;
                cur->y = y * 10.0f;
                y += 1.0f;
                cur = cur->next;
            }
            x += 15.0f;
        }
    } else {
        fprintf(stderr, "place tallest\n");
        /* Place nodes in the tallest row */
        y = 0.0f;
        x = 0.0f;
        cur = nodeRows[tallestRow];
        while(cur != NULL) {
            y = starty;
            x = starty;
            cur->x = x;
            cur->y = y;
            y -= 1.0f;
            cur = cur->next;
        }

        fprintf(stderr, "propagate backwards\n");
        x = -15.0f;
        /* Propagate placement backwards */
        for(int i = tallestRow - 1; i >= 0; i--) {
            cur = nodeRows[i];
            while(cur != NULL) {
                for(auto it = cur->sibs.begin(); it != cur->sibs.end(); it++) {
                    cur->y += (*it)->y;
                    cur->weight ++;
                }
                cur->x = x;
                cur->y /= cur->weight;

                cur = cur->next;
            }
            x -= 15.0f;
        }

        fprintf(stderr, "propagate forwards\n");
        x = 15.0f;
        /* Propagate placement forwards */
        for(int i = tallestRow; i < nodeRowCount - 1; i++) {
            cur = nodeRows[i];
            while(cur != NULL) {
                for(auto it = cur->sibs.begin(); it != cur->sibs.end(); it++) {
                    (*it)->y += cur->y;
                    (*it)->weight ++;
                }
                cur = cur->next;
            }
            cur = nodeRows[i + 1];
            while(cur != NULL) {
                cur->x = x;
                cur->y /= cur->weight;

                cur = cur->next;
            }
            x += 15.0f;
        }
    }
}

void writeNodes()
{
    Node *cur;

    for(int row = 0; row < nodeRowCount; row++) {
	cur = nodeRows[row];
	while(cur != NULL) {
	    if(cur->name == "loss") {
		cur = cur->next;
		continue;
	    }
            for(auto it = cur->sibs.begin(); it != cur->sibs.end(); it++) {
                Node *sib = *it;
		if(sib->name == "loss") {
		    printf("line %f %f %f %f 1 0 0 3\n",
			cur->x, cur->y, cur->x + 5, cur->y + 5);
		    printf("line %f %f %f %f 1 0 0 3\n",
			cur->x + 7, cur->y + 7, cur->x + 9, cur->y + 9);
		} else {
		    printf("line %f %f %f %f 1 1 1 1\n",
			cur->x, cur->y, sib->x, sib->y);
		}
	    }

	    printf("point %f %f 0 1 0 4\n", cur->x, cur->y);
	    printf("string %f %f 1 1 1 %s\n", cur->x, cur->y, cur->name.c_str());
	    cur = cur->next;
	}
    }
}

int main(int argc, char **argv)
{
    readNodes();

    orderNodesInRows();

    placeNodes();

    writeNodes();
}
