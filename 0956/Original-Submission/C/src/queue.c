#include "pampac.h"

/**********************************************************************/
/* Routines implementing a queue in C using linked lists              */
/**********************************************************************/

/*initialize a queue*/
void init_queue (Queue *q) {
  q->head = NULL;
  q->tail = NULL;
}

/*insert an element to the end of the queue*/
void enqueue (Queue *q, PTnode *_value) {
  //allocate a new QueueElement for _value
  QueueElement *newElement;
  newElement = (QueueElement*) malloc (sizeof (QueueElement));
  newElement->value = _value;
  newElement->next = NULL;
  if (q->head == NULL) {
    //first element
    q->head = newElement;
    q->tail = newElement;
  } else {
    //put it to the tail
    q->tail->next = newElement;
    q->tail = newElement;
  }
}

/*delete the first element from the queue*/
void dequeue (Queue *q) {
  QueueElement *element;
  if (q->head == NULL) {
    //empty queue
    return;
  } else {
    element = q->head;
    q->head = q->head->next;
    free (element);
  }
}

/*get the front value of the queue, but don't delete it*/
PTnode* front_of_queue (Queue *q) {
  return q->head->value;
}

/*check if the queue is empty*/
int empty_queue (Queue *q) {
  return (q->head == NULL ? 1:0);
}
