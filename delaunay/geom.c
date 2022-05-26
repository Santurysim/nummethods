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

node_t *find(point_t *pt, point_t *target)
{
    node_t *result;
    if (!pt->head)
        return NULL;
    if (pt->head->next == pt->head) {
        if (pt->head->point == target)
            return pt->head;
        else
            return NULL;
    }
    for (result = pt->head->next; result != pt->head; result = result->next)
        if (result->point == target)
            return result;
    return NULL;
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

int is_right_of(point_t const *pt, segment_t const *s)
{
    double tmp = (s->b->x - s->a->x) * (s->b->y - pt->y)
               - (s->b->y - s->a->y) * (s->b->x - pt->x);
    return (tmp < 0.0);
}

int is_left_of(point_t const *pt, segment_t const *s)
{
    double tmp = (s->b->x - s->a->x) * (s->b->y - pt->y)
               - (s->b->y - s->a->y) * (s->b->x - pt->x);
    return (tmp > 0.0);
}

segment_t lower_common_tangent(node_t *ch1, node_t *ch2)
{
    segment_t s;
    node_t *tmp;
    node_t *x = ch_rightmost_point(ch1);
    node_t *y = ch_leftmost_point(ch2);
    point_t *z = y->next->point; // TODO
    point_t *ztick = x->next->point; // TODO
    point_t *zticktick = find(x, ztick)->previous->point;
    s.a = x->point;
    s.b = y->point;
    for (;;) {
        if (is_right_of(z, &segment)) {
            tmp = find(z, y->point)->next;
            z = y->point;
            y = tmp;
            s.b = y->point;
        } else if (is_right_of(zticktick, &segment)) {
            tmp = fund(zticktick, x->point)->previous;
            zticktick = x->point;
            x = tmp;
            s.a = x->point;
        } else return s;
    }
}

segment_t upper_common_tangent(node_t *ch1, node_t *ch2)
{
    segment_t s;
    node_t *tmp;
    node_t *x = ch_rightmost_point(ch1);
    node_t *y = ch_leftmost_point(ch2);
    point_t *z = y->next->point; // TODO
    point_t *ztick = x->next->point; // TODO
    point_t *zticktick = find(x, ztick)->previous->point;
    s.a = x->point;
    s.b = y->point;
    for (;;) {
        if (is_left_of(z, &segment)) {
            tmp = find(z, y->point)->next;
            z = y->point;
            y = tmp;
            s.b = y->point;
        } else if (is_left_of(zticktick, &segment)) {
            tmp = fund(zticktick, x->point)->previous;
            zticktick = x->point;
            x = tmp;
            s.a = x->point;
        } else return s;
    }
}

int qtest(point_t *h, point_t *i, point_t *j, point_t *k)
{
    double a = h->x * i->y + h->y * j->x + i->x * j->y
             - i->y * j->x - h->y * i->x - h->x * j->y;
    double b = SQRNORM(h) * i->y + h->y * SQRNORM(j) + SQRNORM(i) * j->y
             - i->y * SQRNORM(j) - h->y * SQRNORM(i) - SQRNORM(h) * j->y;
    double c = SQRNORM(h) * i->x + h->x * SQRNORM(j) + SQRNORM(i) * j->x
             - i->x * SQRNORM(j) - h->x * SQRNORM(i) - SQRNORM(h) * j->x;
    double d = SQRNORM(h) * i->x * j->y + h->x * i->y * SQRNORM(j)
             + h->y * SQRNORM(i) * j->x
             - h->y * i->x * SQRNORM(j) - h->x * SQRNORM(i) * j->y
             - SQRNORM(h) * i->y * j->x;
    if (a > 0.0)
        return (a * SQRNORM(k) - b * k->x + c * k->y - d) >= 0.0;
    else
        return (a * SQRNORM(k) - b * k->x + c * k->y - d) <= 0.0;
}

double dist2 (point_t const *a, point_t const *b) {
    return sqrt(SQR(b->x - a->x) + SQR(b->y - a->x));
}

void merge(segment_t ut, segment_t bt)
{
    point_t *r1, *r2, *l1, *l2;
    segment_t s;
    point_t *l = bt.a;
    point_t *r = bt.b;
    while (dist2(ut.a, bt.a) < EPS && dist2(ut.b, bt.b)) {
        int a = 0, b = 0;
        insert_point(l, r);
        insert_point(r, l);
        r1 = find(r, l)->previous;
        s.a = l;
        s.b = r;
        if (is_left_of(r1, &s)) {
            r2 = find(r, r1)->previous;
            while (qtest(r1, l, r, r2)) {
                delete_point(r, r1);
                delete_point(r1, r);
                r1 = r2;
                r2 = find(r, r1)->previous;
            }
        } else
            a = 1;
        l1 = find(l, l1)->next;
        s.a = r;
        s.b = l;
        if (is_right_of(l1, s)) {
            l2 = find(l, l1)->next;
            while (qtest(l, r, l1, l2)) {
                delete_point(l, l1);
                delete_point(l1, l);
                l1 = l2;
                l2 = find(l, l1);
            }
        } else
            b = 1;
        if (a)
            l = l1;
        else
            if (b)
                r = r1;
            else
                if (qtest(l, r, r1, l1))
                    r = r1;
                else
                    l = l1;
        bt.a = l;
        bt.b = r;
    }
}

void trianglulate (point_t *points)
{
    
}
