/* integration/append.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

static inline
void append_interval (gsl_integration_workspace * workspace,
                      double a1, double b1, double area1, double error1)
{
  const size_t i_new = workspace->size ;

  workspace->alist[i_new] = a1;
  workspace->blist[i_new] = b1;
  workspace->rlist[i_new] = area1;
  workspace->elist[i_new] = error1;
  workspace->order[i_new] = i_new;
  workspace->level[i_new] = 0;

  workspace->size++;
}

