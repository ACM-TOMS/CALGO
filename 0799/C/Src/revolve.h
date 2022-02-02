#ifndef _REVOLVE_H_
#define _REVOLVE_H_

enum action { advance, takeshot, restore, firsturn, youturn, terminate, error};

int maxrange(int ss, int tt);
int adjustsize(int* steps, int* snaps, int* reps);
enum action revolve(int* check,int* capo,int* fine,int snaps,int* info);

#endif
