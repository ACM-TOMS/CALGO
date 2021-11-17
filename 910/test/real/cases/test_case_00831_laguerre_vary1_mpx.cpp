
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00831_laguerre_vary1_mpx : public TestCaseReal
    {
    public:
      TestCase_case_00831_laguerre_vary1_mpx() { }
      virtual ~TestCase_case_00831_laguerre_vary1_mpx() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00831_laguerre_vary1_mpx");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const e_float ak = ef::third()   + static_cast<INT32>(100 - (10 * k));
          const e_float bk = ef::quarter() + static_cast<INT32>(0   + (10 * k));
          data[k] = ef::laguerre(ak, bk, ef::euler_gamma() + (100 * k));
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("0.029830924120183453219395462589656282497007108042787847654955022927529664824276444829401147219629611195273029821467559247855171482180625113205893359643685350120336923600247900846358909042752336622672581054223116782488263900836157013028947241099796686421806698987162516866032607701943847241273481587728514696814180658991712662135812116568677475206548962505230561779991557423339137619023135998182215491327"),
           e_float("1.1315826274074035851417337245023754085896834328394593543070178321880575765118676550902542470177371205158001244420600692930913030463701275633639875508936311910292433549814737058671696895292335071710654246587676474008745666468937560590801971102261787623789351065635173647288815347094037420945471592830087010287210459885048887147361287623119417766067415570093447913788945148565743484231222347420623407107e20"),
           e_float("-6.904236466851245676651901851048262215713224630609699366947754799486419527654240809946781585721575976820733583291344666432894214182606207361718032546836199302942199912154634091028683341391718578619272917008319986211359900017235092729535960350456730894424512130529955136470091107147301729788760671877872885599634480489243301474642012853081160138516037145993429962965925144515904500860189912101844250495e38"),
           e_float("-7.7098823816997579929313624670879355265563307391664263304886901328559345991823012708174621739251849583956252604025022869236367996280623164813943106541754303882545987059389300331992341395080819565498229793675368720181863220238785524129864490126818782481415580206018820510543865003227437065104315089714264220444538907952211772310607889833628309534722195461034537198473094582722590385258905210711461446012e55"),
           e_float("-1.5372568830656681964690432239948153114942028783714877005110182340796947713225594639512288160037572702386471184230163492315665828491708671195974557346641706622071783490430167236515519738045132708192638920235761948490383523498564342020778535017963848927264624212571026210153120739038426485154646407750173348078647286149152372131988565692222659809684322331669867923574241794077994699069517139276143758911e77"),
           e_float("-1.7363363333250360647593877302118369405955570228450692372937659194984023443956805875635253793245127652464652979054641504333822616687254690573334320738495249998133089783792337467831354716089287298120946997612725694424051959309565502957121212852994962058458319435843754303968167330995110221923918002349596622039874489038950051087917723118378745658818449013430735134945601534509886647097772645630875494892e107"),
           e_float("-4.4094546568853716385873027507445649043004713765062012097093355041936189621380081129327614741089152815594608446115931687536636932852003131508379572300649841818048748685910552407807138577242266398461044073719231987945021441674969175941720159619939474523459140892708016355469835554010960332300854584369063079166855798569127342016077331083970178023698476200518492631232603976926452502657969793289994695463e140"),
           e_float("-9.5950052733823584105633670502514623116450622264551726167656131806427099567455973780243505956991907841668437904192436267233792256804270249242529156443427090413643294307508325707169018389998854639636408699341621647906847225879238931352176203727299061467203672891581397105109560012152809887300526359770795838349223030262657544012870003240530117077743500602291725837949026619654609643225215326336256677325e175"),
           e_float("-4.0478190404885372304647723111060447470640512897782742659011392414874612273654383878833989460569586929384772802683114461197004708524342862065703082913103761864136201423016062299765178700345621262717325952509448727964180800984695340583299513766162441004827299921610655593688215724945673210604477698375394352682130486041035301440366762927032755647833645166924926419278627319631541619805406339044131917792e212"),
           e_float("-1.4409555055084820363914205918723932420117441395193773899530856639532832912614846685913394081848277374323318205792904650657941074805388200288603904213577009837937582882220091293526499718987567588439083237739745039968924807841859284275927127125976970664046005014758987599720170243847170688755035827747352890377425519061983861107675838550237367050762461226469373127835249837461531420518764593102986028787e250"),
           e_float("-2.5753259739767372234292416873301852996248650739671827386233439190059296208316703413808148019952627608343298081099003700139352756128968232861906457612723566749413161299970179165549825664466022566781553441938974139762893394325552088961490050248992952920738157657417591514146367000241212605051590323329890695561372794611915674479891216107118888249410178378259774792775820158028245288212887002103698974481e288"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00831_laguerre_vary1_mpx(const bool b_write_output)
    {
      return TestCase_case_00831_laguerre_vary1_mpx().execute(b_write_output);
    }
  }
}
