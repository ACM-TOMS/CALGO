with text_io,Solutions;              use text_io,Solutions;
with Complex_Polynomial_Systems;     use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;  use Complex_Polynomial_Systems_io;
with PHCPACK;

procedure use_phc is

  infile,outfile : file_type;        -- input and output file
  p,q : Link_to_Poly_Sys;            -- target and start system
  mixed_volume : natural;            -- root count is mixed volume
  sols : Solution_List;              -- list of solutions

begin
  Open(infile,in_file,"test.in");
  get(infile,p);
  Create(outfile,out_file,"test.out");
  put(outfile,p.all);
  q := new Poly_Sys(p'range);
  PHCPACK.Static_Lifting(outfile,p.all,mixed_volume,q.all,sols);
  PHCPACK.Artificial_Parameter_Continuation(outfile,p.all,q.all,sols);
  PHCPACK.Refine_Roots(outfile,p.all,sols);
end use_phc;
