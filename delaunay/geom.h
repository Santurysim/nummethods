#ifndef GEOM_H
#define GEOM_H

#define EPS 1e-16

struct point;
struct node;
struct adj_list;
struct segment;

typedef struct point point_t;
typedef struct node node_t;
typedef struct adj_list adj_list_t;
typedef struct segment segment_t;

struct point {
    double x;
    double y;
    node_t *head;
};

struct node {   // Weird double linked list
    point_t *point;
    node_t *next;
    node_t *previous;
};

struct segment {
    point_t *a;
    point_t *b;
};

void insert_point(point_t *pt, point_t *newpt);

void delete_point(point_t *pt, point_t *oldpt);

node_t *ch_rightmost_point(node_t *ch);
node_t *ch_leftmost_point(node_t *ch);

segment_t lowest_common_tangent(node_t *ch1, node_t *ch2);

int qtest(point_t *, point_t *, point_t *, point_t *);

void merge(/* TODO */);

#endif /* GEOM_H */
