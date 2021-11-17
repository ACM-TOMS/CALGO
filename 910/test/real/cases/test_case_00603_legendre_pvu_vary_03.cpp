
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00603_legendre_pvu_vary_03 : public TestCaseReal
    {
    public:
      TestCase_case_00603_legendre_pvu_vary_03() { }
      virtual ~TestCase_case_00603_legendre_pvu_vary_03() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00603_legendre_pvu_vary_03");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        static const e_float delta = ef::euler_gamma() / 100;
        static const e_float sqrt_1_3 = ef::sqrt(ef::third());
        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const e_float xk = ef::one_minus() +  (e_float(k) / 5);
          const e_float x  = k <= 5 ? xk + delta : xk - delta;
          const e_float v  = +((10 - k) * 13) + sqrt_1_3;
          const e_float u  = -((10 - k) * 11) - sqrt_1_5;
          data[k] = ef::legendre_p(v, u, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("-1.1484194954389099564013088332505347484137095748318446488388276292828381314895871766392523975374543961675901142642517835754212841704201516064827561261854071944099164745048210481221180259867118647610929250534466691330525296901356794733957112099110419631838332447406552040610223668167914694353507909369133538654289551434539886194441875464357236785382516184999543051285331922156631175926108562375280629382e-136"),
           e_float("-1.5133046201576355179826841297249706029298011246953138688346052323419986393300539754706479691427159394526722067433376254556437689130523681134790741527976073066221988802189864829324592430753998498263126434021858400015071757594117114092934758234312740517299390336743741085330038711251922416411568512725602749839108424537104079913386944963091515068231546620427578114233671898963599640730752150160918149328e-191"),
           e_float("-4.9499342942779953239021081377144345568107497148348726307628475376788446643289725614610546496447497020430992351142813925995708611633918571290013951695506375737085099297587111267855368966884594697210189552280478256625159133961947178142265752785196981297426445742508536444843266324066462068194675233085909202047475332425039492506686901169105831030433425465881510270599349242177553520689959649211625433518e-174"),
           e_float("-2.1785494204289742844883482618377793560618698535158447131331467630481803025168697104710416453912332880524792993944023764913599368425704534982978388192578210734518075938465390507949852400326716644720030988160568822769730199177423751333700393717416195132446303601482291403351972400619497031572888627286646641295100135961516953401753679871870642109787830373512532624377689780710136212551628631187092524057e-148"),
           e_float("-1.9328579733417161031321039815721154924142656120729126904529327188003968465761618344764568131660688005649163989855519172384426100543273792722157176826727458352603631835864895158637888886533835533029998399184589085123980423235694848215297441742656330334312930049437099784335991719465441715456009674941162372791554956844096026016566065349388725227649655460863264516092163768728557291834077143186224230447e-123"),
           e_float("-9.7915111811389520049795792533200166383831193846946138585597786687240915049823822718332953865683366997627951235727823332081751945287448407740327009777955188632816052982381291823665469601568345239234758608058902082916598398670717795823484470795040285948767255216320767787965400817397645422332503040605201818487933356897698553476991716859963401535621089946732496414984284798598314386009080723320187945888e-99"),
           e_float("2.1222734536408388029983703455651754681498219552504058653756385592639008596990232737782733921772018246445347701675416078084373440716470812313149581582498283393230813921496895646765519897080675823726508371895668461696537316616672730961355287814295120903799031442546634336509208035148078203452496465727053849888540068152465478047184240540852474631734294159443032900220118947600220630182345608711659294315e-75"),
           e_float("2.3233322737893137532444429572936992474528073821327924078490404489502679676777468804304114493374493675483255626550675574409317850465069580445546165856365375032024599211795900227421342810432770271868291782575122513212213312885787547078359697479703373963020583531374575547551779706075202415835698697580140759360129890939569530667268947618869616442611984181885887045278681275653303624472656654191152016675e-53"),
           e_float("2.543614544258075362717446101524708203470617059703731870526415699606252963023346520529119231421127147847974147726437560852342023228517742716781656838845916011876553667511353658271825497045444281068937722086924914344772151797042645711374402950294962318337769120581080433481157154382712062449361641801697360013645843151857970690940564148180631931299883450772463016345821028351465035599607370105590918104e-32"),
           e_float("5.9200915661039084952474287900957638342144929043490398562434852472342142419592871484367974585670767279243575661374925061679047364382927961973799719983771994819429184251155604064043495038962425363101400195258846425492273613783198740317282783683939076335843332763540158882859777259450502154720886243999375023797200289657223360009016344269697474070452858862860255437167299181338217571009627162328853417943e-15"),
           e_float("0.30501085812125178033171444746444029465508976414702416117593951153210043046032388821605988696996985227601917011847958606032664794677266570428127505471281149385083365063575160274353653868070316469864716219806992467966120569577992741563864386821568377800673547223379506315409481786860150555785345767207858285323050717571778929132244282605913506574164075804397744925590563384917521141774292214839656479747"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00603_legendre_pvu_vary_03(const bool b_write_output)
    {
      return TestCase_case_00603_legendre_pvu_vary_03().execute(b_write_output);
    }
  }
}
