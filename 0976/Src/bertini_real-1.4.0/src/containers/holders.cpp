#include "containers/holders.hpp"







int PointHolder::add_point(vec_mp new_point)
{
	
	if (num_pts_!=0 && this->pts_mp_==NULL) {
		printf("trying to add point to PointHolder with non-zero num_points and NULL container!\n");
		br_exit(9713);
	}
	
	if (num_pts_==0 && this->pts_mp_!=NULL) {
		printf("trying to add point to PointHolder with num_points==0 and non-NULL container!\n");
		br_exit(9713);
	}
	
	
	if (num_pts_==0) {
		pts_mp_ = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		pts_mp_ = (vec_mp *)br_realloc(pts_mp_, (num_pts_+1) * sizeof(vec_mp));
	}
	
	init_vec_mp2(pts_mp_[num_pts_], new_point->size, new_point->curr_prec);
	pts_mp_[num_pts_]->size = new_point->size;
	vec_cp_mp(pts_mp_[num_pts_], new_point);
	
	num_pts_++;
	
	return num_pts_-1;
}

int PatchHolder::add_patch(vec_mp new_patch)
{
	
	if (num_patches_!=0 && patch_mp_==NULL) {
		printf("trying to add patch to witness set with non-zero num_patches and NULL container!\n");
		deliberate_segfault();
	}
	
	if (num_patches_==0 && patch_mp_!=NULL) {
		printf("trying to add point to witness set with num_points==0 and non-NULL container!\n");
		deliberate_segfault();
	}
	
	
	if (num_patches_==0) {
		patch_mp_ = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		patch_mp_ = (vec_mp *)br_realloc(patch_mp_, (num_patches_+1) * sizeof(vec_mp));
	}
	
	init_vec_mp2(patch_mp_[num_patches_], new_patch->size,new_patch->curr_prec);
	patch_mp_[num_patches_]->size = new_patch->size;
	vec_cp_mp(patch_mp_[num_patches_], new_patch);
	
	this->num_patches_++;
	
	
	int burnsome = 0;
	for (unsigned int ii=0; ii<this->num_patches_; ii++) {
		burnsome += this->patch_mp_[ii]->size;
	}
	

	return num_patches_-1;
}


int LinearHolder::add_linear(vec_mp new_linear)
{
	
	if (num_linears_!=0 && L_mp_==NULL) {
		printf("trying to add linear to linear holder with non-zero num_linears and NULL container!\n");
		br_exit(9711);
	}
	
	if (num_linears_==0 && L_mp_!=NULL) {
		printf("trying to add linear to linear holder with num_linears==0 and non-NULL container!\n");
		br_exit(9711);
	}
	
	
	if (num_linears_==0) {
		L_mp_ = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		L_mp_ = (vec_mp *)br_realloc(L_mp_, (num_linears_+1) * sizeof(vec_mp));
	}
	
	init_vec_mp2(L_mp_[num_linears_], new_linear->size,new_linear->curr_prec);
	L_mp_[num_linears_]->size = new_linear->size;
	vec_cp_mp(L_mp_[num_linears_], new_linear);
	
	this->num_linears_++;
	
	
	return num_linears_-1;
}

