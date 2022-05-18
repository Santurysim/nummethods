#include <stdlib.h>
#include <math.h>

#include "geom.h"

// TODO make these operations work with generic point types
void insert_point(point_t *pt, point_t *newpt)
{
    double newpt_angle;
    node_t *newnode = (node_t*)malloc(sizeof(node_t));
    if (!newnode) exit(1);
    
    newpt_angle = atan2(newpt->y - pt->y, newpt->x - pt->x);
    if (newpt_angle < 0.0) newpt_angle += 2 * M_PI;

    newnode->point = newpt;

    if (!pt->head) {
        pt->head = newnode;
        newnode->next = newnode->previous = newnode;
    } else if (pt->head->next == pt->head) {
        pt->head->next = pt->head->previous = newnode;
        newnode->next = newnode->previous = pt->head;
    } else {
        for (node_t *it = pt->head->next; it != pt->head; it = it->next) {
            double curpt_angle = atan2(it->y - pt->y, it->x - pt->x);
            if (curpt_angle < 0.0) curpt_angle += 2 * M_PI;
            if (curpt_angle < newpt_angle)
                break;
        }
        newnode->previous = it->previous;
        newnode->next = it;
        it->previous->next = newnode;
        it->previous = newnode;
    }
}

void delete_point(point_t *pt, point_t *oldpt)
{
    node_t *oldnode;
    if (pt->head->next == pt->head && pt->head->point == oldpt) {
        free(pt->head);
        pt->head = 0;
        return;
    }
    for (node_t *it = pt->head->next; it != head; it = it->next)
        if (it->point == oldpt) {
            oldnode = it;
            break;
        }
   
    if (oldnode == pt->head)
        pt->head = pt->head->next;
    oldnode->previous->next = oldnode->next;
    oldnode->next->previous = oldnode->previous;
    free(oldnode);
}

node_t *ch_rightmost_point(node_t *ch)
{
    // TODO: why node->next yields right point in a convex hull?
    // Answer: bacause
    node_t *result = ch;
    double xmin = ch->point->x;
    for (node_t *node = ch->next, node != ch; node = node->next) {
        if (node->point->x < xmin) {
            result = node;
            xmin = node->point->x;
        }
    }
    return result;
}

node_t *ch_leftmost_point(node_t *ch)
{
    node_t *result = ch;
    double xmax = ch->point->x;
    for (node_t *node = ch->next, node != ch; node = node->next) {
        if (node->point->x < xmax) {
            result = node;
            xmax = node->point->x;
        }
    }
    return result;
}

segment_t lowest_common_tangent(node_t *ch1, node_t *ch2)
{
    node_t *x = ch_rightmost_point(ch1);
    node_t *y = ch_leftmost_point(ch2);
    point_t *z = x->point->head->point; // TODO
    point_t *ztick = y->point->head->next; // TODO
    for (;;) {

    }
}
