<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile>
  <compound kind="file">
    <name>0_installing.md</name>
    <path>/Users/ongbw/from-vcs/github/ridc/doc/source/</path>
    <filename>0__installing_8md</filename>
  </compound>
  <compound kind="file">
    <name>1_contributing.md</name>
    <path>/Users/ongbw/from-vcs/github/ridc/doc/source/</path>
    <filename>1__contributing_8md</filename>
  </compound>
  <compound kind="file">
    <name>2_running.md</name>
    <path>/Users/ongbw/from-vcs/github/ridc/doc/source/</path>
    <filename>2__running_8md</filename>
  </compound>
  <compound kind="file">
    <name>3_use.md</name>
    <path>/Users/ongbw/from-vcs/github/ridc/doc/source/</path>
    <filename>3__use_8md</filename>
  </compound>
  <compound kind="file">
    <name>brusselator.cpp</name>
    <path>/Users/ongbw/from-vcs/github/ridc/examples/brusselator_gsl/</path>
    <filename>brusselator__gsl_2brusselator_8cpp</filename>
    <includes id="ridc_8h" name="ridc.h" local="yes" imported="no">ridc.h</includes>
    <includes id="brusselator__gsl_2brusselator_8h" name="brusselator.h" local="yes" imported="no">brusselator.h</includes>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>brusselator__gsl_2brusselator_8cpp.html</anchorfile>
      <anchor>a0ddf1224851353fc92bfbff6f499fa97</anchor>
      <arglist>(int argc, char *argv[])</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>brusselator.cpp</name>
    <path>/Users/ongbw/from-vcs/github/ridc/examples/brusselator_mkl/</path>
    <filename>brusselator__mkl_2brusselator_8cpp</filename>
    <includes id="ridc_8h" name="ridc.h" local="yes" imported="no">ridc.h</includes>
    <includes id="brusselator__mkl_2brusselator_8h" name="brusselator.h" local="yes" imported="no">brusselator.h</includes>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>brusselator__mkl_2brusselator_8cpp.html</anchorfile>
      <anchor>a0ddf1224851353fc92bfbff6f499fa97</anchor>
      <arglist>(int argc, char *argv[])</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>brusselator.cpp</name>
    <path>/Users/ongbw/from-vcs/github/ridc/examples/brusselator_radau_mkl/</path>
    <filename>brusselator__radau__mkl_2brusselator_8cpp</filename>
    <includes id="ode_8h" name="ode.h" local="yes" imported="no">ode.h</includes>
    <member kind="function">
      <type>void</type>
      <name>rhs</name>
      <anchorfile>brusselator__radau__mkl_2brusselator_8cpp.html</anchorfile>
      <anchor>a85cca37916ac01b8648f8a864bafdae5</anchor>
      <arglist>(double t, double *u, PARAMETER param, double *f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>newt</name>
      <anchorfile>brusselator__radau__mkl_2brusselator_8cpp.html</anchorfile>
      <anchor>aa04e20511f2366a5ec0bdff1064203d2</anchor>
      <arglist>(double t, double *uprev, double *Kguess, double *g, PARAMETER param, BUTCHER rk)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>jac</name>
      <anchorfile>brusselator__radau__mkl_2brusselator_8cpp.html</anchorfile>
      <anchor>a2f4e36896b92e6193d45ef27b7020780</anchor>
      <arglist>(double t, double *uprev, double *Kguess, double *J, PARAMETER param, BUTCHER rk)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>step</name>
      <anchorfile>brusselator__radau__mkl_2brusselator_8cpp.html</anchorfile>
      <anchor>add04c1b0f2b2529551dd6fa2ec1a2296</anchor>
      <arglist>(double t, double *uold, PARAMETER param, double *unew, BUTCHER rk)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>brusselator.h</name>
    <path>/Users/ongbw/from-vcs/github/ridc/examples/brusselator_gsl/</path>
    <filename>brusselator__gsl_2brusselator_8h</filename>
    <includes id="ridc_8h" name="ridc.h" local="yes" imported="no">ridc.h</includes>
    <class kind="class">Brusselator_GSL</class>
  </compound>
  <compound kind="file">
    <name>brusselator.h</name>
    <path>/Users/ongbw/from-vcs/github/ridc/examples/brusselator_mkl/</path>
    <filename>brusselator__mkl_2brusselator_8h</filename>
    <includes id="ridc_8h" name="ridc.h" local="yes" imported="no">ridc.h</includes>
    <class kind="class">Brusselator_MKL</class>
    <member kind="define">
      <type>#define</type>
      <name>_BRUSSELATOR_H_</name>
      <anchorfile>brusselator__mkl_2brusselator_8h.html</anchorfile>
      <anchor>a812147959eafc576a0f4ebed8ae4d0ca</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>ode.h</name>
    <path>/Users/ongbw/from-vcs/github/ridc/examples/brusselator_radau_mkl/</path>
    <filename>ode_8h</filename>
    <class kind="struct">PARAMETER</class>
    <class kind="struct">BUTCHER</class>
    <member kind="function">
      <type>void</type>
      <name>rhs</name>
      <anchorfile>ode_8h.html</anchorfile>
      <anchor>a85cca37916ac01b8648f8a864bafdae5</anchor>
      <arglist>(double t, double *u, PARAMETER param, double *f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>newt</name>
      <anchorfile>ode_8h.html</anchorfile>
      <anchor>aa04e20511f2366a5ec0bdff1064203d2</anchor>
      <arglist>(double t, double *uprev, double *Kguess, double *g, PARAMETER param, BUTCHER rk)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>jac</name>
      <anchorfile>ode_8h.html</anchorfile>
      <anchor>a2f4e36896b92e6193d45ef27b7020780</anchor>
      <arglist>(double t, double *uprev, double *Kguess, double *J, PARAMETER param, BUTCHER rk)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>step</name>
      <anchorfile>ode_8h.html</anchorfile>
      <anchor>a02a2cfd5bf4a4c35f861e2c2d2dab18e</anchor>
      <arglist>(double t, double *u, PARAMETER param, double *unew, BUTCHER rk)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>radau.cpp</name>
    <path>/Users/ongbw/from-vcs/github/ridc/examples/brusselator_radau_mkl/</path>
    <filename>radau_8cpp</filename>
    <includes id="ode_8h" name="ode.h" local="yes" imported="no">ode.h</includes>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>radau_8cpp.html</anchorfile>
      <anchor>a0ddf1224851353fc92bfbff6f499fa97</anchor>
      <arglist>(int argc, char *argv[])</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>explicit.cpp</name>
    <path>/Users/ongbw/from-vcs/github/ridc/examples/explicit/</path>
    <filename>explicit_8cpp</filename>
    <includes id="ridc_8h" name="ridc.h" local="yes" imported="no">ridc.h</includes>
    <class kind="class">ExplicitOde</class>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>explicit_8cpp.html</anchorfile>
      <anchor>a0ddf1224851353fc92bfbff6f499fa97</anchor>
      <arglist>(int argc, char *argv[])</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>implicit.cpp</name>
    <path>/Users/ongbw/from-vcs/github/ridc/examples/implicit/</path>
    <filename>implicit_8cpp</filename>
    <includes id="ridc_8h" name="ridc.h" local="yes" imported="no">ridc.h</includes>
    <class kind="class">ImplicitOde</class>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>implicit_8cpp.html</anchorfile>
      <anchor>a0ddf1224851353fc92bfbff6f499fa97</anchor>
      <arglist>(int argc, char *argv[])</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>implicit.cpp</name>
    <path>/Users/ongbw/from-vcs/github/ridc/examples/implicit_mkl/</path>
    <filename>mkl_2implicit_8cpp</filename>
    <includes id="ridc_8h" name="ridc.h" local="yes" imported="no">ridc.h</includes>
    <includes id="implicit_8h" name="implicit.h" local="yes" imported="no">implicit.h</includes>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>mkl_2implicit_8cpp.html</anchorfile>
      <anchor>a0ddf1224851353fc92bfbff6f499fa97</anchor>
      <arglist>(int argc, char *argv[])</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>implicit.h</name>
    <path>/Users/ongbw/from-vcs/github/ridc/examples/implicit_mkl/</path>
    <filename>implicit_8h</filename>
    <includes id="ridc_8h" name="ridc.h" local="yes" imported="no">ridc.h</includes>
    <class kind="class">ImplicitMKL</class>
    <member kind="define">
      <type>#define</type>
      <name>_IMPLICIT_H_</name>
      <anchorfile>implicit_8h.html</anchorfile>
      <anchor>ae86b67c4f0c0cc9a5f89ed888e72b8d0</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>README.md</name>
    <path>/Users/ongbw/from-vcs/github/ridc/</path>
    <filename>README_8md</filename>
  </compound>
  <compound kind="file">
    <name>ridc.cpp</name>
    <path>/Users/ongbw/from-vcs/github/ridc/src/</path>
    <filename>ridc_8cpp</filename>
    <includes id="ridc_8h" name="ridc.h" local="yes" imported="no">ridc.h</includes>
    <member kind="function">
      <type>void</type>
      <name>ridc_fe</name>
      <anchorfile>ridc_8cpp.html</anchorfile>
      <anchor>a088c5b6cd62faca40e47d5a428a0877a</anchor>
      <arglist>(ODE *ode, int order, double *sol)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ridc_be</name>
      <anchorfile>ridc_8cpp.html</anchorfile>
      <anchor>a44576f55de3e48c9752f1e75f31b1c72</anchor>
      <arglist>(ODE *ode, int order, double *sol)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>lagrange_coeff</name>
      <anchorfile>ridc_8cpp.html</anchorfile>
      <anchor>a73647697e05b092fe6d4a43361977b0a</anchor>
      <arglist>(double *x, int Nx, int i, double *L)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_quad_weight</name>
      <anchorfile>ridc_8cpp.html</anchorfile>
      <anchor>a07e5dce091c0ec540e8f9aff1bbc2a66</anchor>
      <arglist>(double *L, int Nx, double a, double b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>integration_matrices</name>
      <anchorfile>ridc_8cpp.html</anchorfile>
      <anchor>a4bd004e1c1236beb766757b596fae927</anchor>
      <arglist>(int N, double **S)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_unif_nodes</name>
      <anchorfile>ridc_8cpp.html</anchorfile>
      <anchor>a6784be1bca782a9aef5584921af9304c</anchor>
      <arglist>(double *x, int Nx, double a, double b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>corr_fe</name>
      <anchorfile>ridc_8cpp.html</anchorfile>
      <anchor>a9d505d3e3fe2d373768b906cc207b199</anchor>
      <arglist>(ODE *ode, double *uold, double **fprev, double **S, int index, int level, double t, double *unew)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>corr_be</name>
      <anchorfile>ridc_8cpp.html</anchorfile>
      <anchor>a8f0874733447d234c1f21129a19c5630</anchor>
      <arglist>(ODE *ode, double *uold, double **fprev, double **S, int index, int level, double t, double *unew)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>ridc.h</name>
    <path>/Users/ongbw/from-vcs/github/ridc/src/</path>
    <filename>ridc_8h</filename>
    <class kind="class">ODE</class>
    <member kind="function">
      <type>void</type>
      <name>ridc_fe</name>
      <anchorfile>ridc_8h.html</anchorfile>
      <anchor>a088c5b6cd62faca40e47d5a428a0877a</anchor>
      <arglist>(ODE *ode, int order, double *sol)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ridc_be</name>
      <anchorfile>ridc_8h.html</anchorfile>
      <anchor>a44576f55de3e48c9752f1e75f31b1c72</anchor>
      <arglist>(ODE *ode, int order, double *sol)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>lagrange_coeff</name>
      <anchorfile>ridc_8h.html</anchorfile>
      <anchor>a73647697e05b092fe6d4a43361977b0a</anchor>
      <arglist>(double *x, int Nx, int i, double *L)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_quad_weight</name>
      <anchorfile>ridc_8h.html</anchorfile>
      <anchor>a07e5dce091c0ec540e8f9aff1bbc2a66</anchor>
      <arglist>(double *L, int Nx, double a, double b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>integration_matrices</name>
      <anchorfile>ridc_8h.html</anchorfile>
      <anchor>a6ae6b769bd00841b8206f3c73b7e9ec2</anchor>
      <arglist>(int Nx, double **S)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init_unif_nodes</name>
      <anchorfile>ridc_8h.html</anchorfile>
      <anchor>a6784be1bca782a9aef5584921af9304c</anchor>
      <arglist>(double *x, int Nx, double a, double b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>corr_fe</name>
      <anchorfile>ridc_8h.html</anchorfile>
      <anchor>a9d505d3e3fe2d373768b906cc207b199</anchor>
      <arglist>(ODE *ode, double *uold, double **fprev, double **S, int index, int level, double t, double *unew)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>corr_be</name>
      <anchorfile>ridc_8h.html</anchorfile>
      <anchor>a8f0874733447d234c1f21129a19c5630</anchor>
      <arglist>(ODE *ode, double *uold, double **fprev, double **S, int index, int level, double t, double *unew)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Brusselator_GSL</name>
    <filename>classBrusselator__GSL.html</filename>
    <base>ODE</base>
    <member kind="function">
      <type></type>
      <name>Brusselator_GSL</name>
      <anchorfile>classBrusselator__GSL.html</anchorfile>
      <anchor>a8479c3c3837aa34e6d8a44e40305b0e6</anchor>
      <arglist>(int my_neq, int my_nt, double my_ti, double my_tf, double my_dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rhs</name>
      <anchorfile>classBrusselator__GSL.html</anchorfile>
      <anchor>abee5887ab8e67d01da6bc4ff646dd572</anchor>
      <arglist>(double t, double *u, double *f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>step</name>
      <anchorfile>classBrusselator__GSL.html</anchorfile>
      <anchor>a8ed20dea95f4b1705030f0aeee2deacb</anchor>
      <arglist>(double t, double *u, double *unew)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>neq</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>ad10440423b2185d322223da17d1135c5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nt</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a275faaf08a8602b6a0c60131d3b874b0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ti</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>abaa6d4b370ec903c7a1f53de5d758cf1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>tf</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>aa53beae63fb47d22abca58dd4a407e9e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dt</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a489591849cd00a583407cde072b51acd</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>newt</name>
      <anchorfile>classBrusselator__GSL.html</anchorfile>
      <anchor>adfbabbb536eec4c8c7c8ea0e6e04b84f</anchor>
      <arglist>(double t, double *uprev, double *uguess, double *g)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>jac</name>
      <anchorfile>classBrusselator__GSL.html</anchorfile>
      <anchor>ae44dc7d932b081c1a06359aab206dfb7</anchor>
      <arglist>(double t, double *u, double *J)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Brusselator_MKL</name>
    <filename>classBrusselator__MKL.html</filename>
    <base>ODE</base>
    <member kind="function">
      <type></type>
      <name>Brusselator_MKL</name>
      <anchorfile>classBrusselator__MKL.html</anchorfile>
      <anchor>a356b6fb1a3a2a6f9e79795b1834b3c71</anchor>
      <arglist>(int my_neq, int my_nt, double my_ti, double my_tf, double my_dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rhs</name>
      <anchorfile>classBrusselator__MKL.html</anchorfile>
      <anchor>a90e55526052241ef810c61cf7ebf382f</anchor>
      <arglist>(double t, double *u, double *f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>step</name>
      <anchorfile>classBrusselator__MKL.html</anchorfile>
      <anchor>a3193b7fb7ec213ebd468eafa1e66c020</anchor>
      <arglist>(double t, double *u, double *unew)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>neq</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>ad10440423b2185d322223da17d1135c5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nt</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a275faaf08a8602b6a0c60131d3b874b0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ti</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>abaa6d4b370ec903c7a1f53de5d758cf1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>tf</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>aa53beae63fb47d22abca58dd4a407e9e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dt</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a489591849cd00a583407cde072b51acd</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>newt</name>
      <anchorfile>classBrusselator__MKL.html</anchorfile>
      <anchor>ab4859e93248da4c28f7d06aa4562675c</anchor>
      <arglist>(double t, double *uprev, double *uguess, double *g)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>jac</name>
      <anchorfile>classBrusselator__MKL.html</anchorfile>
      <anchor>ace16b06feefec053278637a23cb94bad</anchor>
      <arglist>(double t, double *u, double *J)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>BUTCHER</name>
    <filename>structBUTCHER.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>S</name>
      <anchorfile>structBUTCHER.html</anchorfile>
      <anchor>a175f670f536987c25dcfdcb0856b1393</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>b</name>
      <anchorfile>structBUTCHER.html</anchorfile>
      <anchor>a294a17cbbdd7e1f33f2fcb7d912c1bb9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double *</type>
      <name>c</name>
      <anchorfile>structBUTCHER.html</anchorfile>
      <anchor>a3dc46ec7591b13f9ba6d948e10729d62</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double **</type>
      <name>A</name>
      <anchorfile>structBUTCHER.html</anchorfile>
      <anchor>a4b7de674ab073c62428dffd65f15a0d2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>ExplicitOde</name>
    <filename>classExplicitOde.html</filename>
    <base>ODE</base>
    <member kind="function">
      <type></type>
      <name>ExplicitOde</name>
      <anchorfile>classExplicitOde.html</anchorfile>
      <anchor>a706cb79b8b7fc8518d985f670664d1b3</anchor>
      <arglist>(int my_neq, int my_nt, double my_ti, double my_tf, double my_dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rhs</name>
      <anchorfile>classExplicitOde.html</anchorfile>
      <anchor>a6298c3cda439d026165d08d1c1858ab4</anchor>
      <arglist>(double t, double *u, double *f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>step</name>
      <anchorfile>classExplicitOde.html</anchorfile>
      <anchor>af0d28f83df1cbef1a92433389341597b</anchor>
      <arglist>(double t, double *u, double *unew)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>neq</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>ad10440423b2185d322223da17d1135c5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nt</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a275faaf08a8602b6a0c60131d3b874b0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ti</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>abaa6d4b370ec903c7a1f53de5d758cf1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>tf</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>aa53beae63fb47d22abca58dd4a407e9e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dt</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a489591849cd00a583407cde072b51acd</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>ImplicitMKL</name>
    <filename>classImplicitMKL.html</filename>
    <base>ODE</base>
    <member kind="function">
      <type></type>
      <name>ImplicitMKL</name>
      <anchorfile>classImplicitMKL.html</anchorfile>
      <anchor>a787875b2e52713705ba39c56f68f072b</anchor>
      <arglist>(int my_neq, int my_nt, double my_ti, double my_tf, double my_dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rhs</name>
      <anchorfile>classImplicitMKL.html</anchorfile>
      <anchor>ad5d1c89afc562a0185346ff8e1fe909f</anchor>
      <arglist>(double t, double *u, double *f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>step</name>
      <anchorfile>classImplicitMKL.html</anchorfile>
      <anchor>a69936580278402c177332149d044d016</anchor>
      <arglist>(double t, double *u, double *unew)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>neq</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>ad10440423b2185d322223da17d1135c5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nt</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a275faaf08a8602b6a0c60131d3b874b0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ti</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>abaa6d4b370ec903c7a1f53de5d758cf1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>tf</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>aa53beae63fb47d22abca58dd4a407e9e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dt</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a489591849cd00a583407cde072b51acd</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>newt</name>
      <anchorfile>classImplicitMKL.html</anchorfile>
      <anchor>ac9153700f10f89e31b0db45bd2c56c66</anchor>
      <arglist>(double t, double *uprev, double *uguess, double *g)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>jac</name>
      <anchorfile>classImplicitMKL.html</anchorfile>
      <anchor>aea67cdda27bfd080c4914f92bfa566c7</anchor>
      <arglist>(double t, double *u, double *J)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>ImplicitOde</name>
    <filename>classImplicitOde.html</filename>
    <base>ODE</base>
    <member kind="function">
      <type></type>
      <name>ImplicitOde</name>
      <anchorfile>classImplicitOde.html</anchorfile>
      <anchor>a39bc33f4917328083edcf86d1f61c7d5</anchor>
      <arglist>(int my_neq, int my_nt, double my_ti, double my_tf, double my_dt)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rhs</name>
      <anchorfile>classImplicitOde.html</anchorfile>
      <anchor>a19395eeb6ddb02ef64621b5019bbaf53</anchor>
      <arglist>(double t, double *u, double *f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>step</name>
      <anchorfile>classImplicitOde.html</anchorfile>
      <anchor>a903918de36451c9bf72b0bce1561d5b2</anchor>
      <arglist>(double t, double *u, double *unew)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>neq</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>ad10440423b2185d322223da17d1135c5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nt</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a275faaf08a8602b6a0c60131d3b874b0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ti</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>abaa6d4b370ec903c7a1f53de5d758cf1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>tf</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>aa53beae63fb47d22abca58dd4a407e9e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dt</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a489591849cd00a583407cde072b51acd</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>ODE</name>
    <filename>classODE.html</filename>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>rhs</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a9499def749f2914a41b4616f879c991b</anchor>
      <arglist>(double t, double *u, double *f)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>step</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a966f35008ac30511d950b557f53b7468</anchor>
      <arglist>(double t, double *u, double *unew)=0</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>neq</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>ad10440423b2185d322223da17d1135c5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nt</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a275faaf08a8602b6a0c60131d3b874b0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ti</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>abaa6d4b370ec903c7a1f53de5d758cf1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>tf</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>aa53beae63fb47d22abca58dd4a407e9e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dt</name>
      <anchorfile>classODE.html</anchorfile>
      <anchor>a489591849cd00a583407cde072b51acd</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>PARAMETER</name>
    <filename>structPARAMETER.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>neq</name>
      <anchorfile>structPARAMETER.html</anchorfile>
      <anchor>a8ee37a1980633da35df1c8abfb1ec749</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nt</name>
      <anchorfile>structPARAMETER.html</anchorfile>
      <anchor>a2d5f4361422666a88a03950206b5a3ff</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ti</name>
      <anchorfile>structPARAMETER.html</anchorfile>
      <anchor>af5efb24886ce50393018b9ccab35f45e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>tf</name>
      <anchorfile>structPARAMETER.html</anchorfile>
      <anchor>a59b3b90e4ce6d567e0abd94e52f8aa1b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dt</name>
      <anchorfile>structPARAMETER.html</anchorfile>
      <anchor>a6c8d30eb1e07202b09f07c946753e482</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>md_README</name>
    <title>About RIDC</title>
    <filename>md_README</filename>
  </compound>
  <compound kind="page">
    <name>md_doc_source_0_installing</name>
    <title>Building and Installing</title>
    <filename>md_doc_source_0_installing</filename>
  </compound>
  <compound kind="page">
    <name>page_contributing</name>
    <title>Contributing</title>
    <filename>page_contributing</filename>
  </compound>
  <compound kind="page">
    <name>running</name>
    <title>Running the Examples</title>
    <filename>running</filename>
  </compound>
  <compound kind="page">
    <name>md_doc_source_3_use</name>
    <title>Using the RIDC Library</title>
    <filename>md_doc_source_3_use</filename>
  </compound>
</tagfile>
