
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00685_legendre_pnu_vary_05 : public TestCaseReal
    {
    public:
      TestCase_case_00685_legendre_pnu_vary_05() { }
      virtual ~TestCase_case_00685_legendre_pnu_vary_05() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00685_legendre_pnu_vary_05");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        static const e_float delta = ef::euler_gamma() / 100;
        static const e_float sqrt_1_5 = ef::sqrt(ef::fifth());
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const e_float xk = ef::one_minus() +  (e_float(k) / 5);
          const e_float x  = k <= 5 ? xk + delta : xk - delta;
          const INT32   n  = +127;
          const e_float u  = +127 + sqrt_1_5;
          data[k] = ef::legendre_p(n, u, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("-4.4684896012599000021174466847307963907990044425727857043347673168610062472203577051898042203667480524516973369180704611707069301494025503943646651680683452572756487643316342059876096363626536111677821861184310843454833791164586440538110359165148083396773856487153626120434909954901295889113797258201455584360800781260144012976160265949018876449950889934245310872552741162688794210825735975157250590597e126"),
           e_float("-4.6866680154817777968906207843595800589933987825597109630600095089646816247007895267180051483396354596158371529337093476463989670914298417543741068272349655285660980283388057140771463894523687582749918439450441188597006362544443782367140183635052209069987339489544786064035049677756065775037489263844879000515657838445433860943130926421970576322397079559947093343691877111752013710240664805588060042371e222"),
           e_float("-1.7684734278076715156542231189376400178349906015810102506984131270003441177054428284503437255369378107617467828338643616172194797325987102758991522028321107203012537435514403257994026547670731210938611564887837654405817462716075707361477179919482798520703305803459310838834264278454756131364998495405579970917899772341075002773025998420354944849925815699271855235683814361357760609639272452254753634939e238"),
           e_float("-5.0731857293906888602975605609187746627982100843506862486264165680773431737210898891002132702041673690130229058189042249332180770894942267531848065548512345458346427059847147505414830084972109949220877995903570759236854308132675608260287500578034941153205844763156066736082289037604804308503970420642820410982584078807256647413108493500453832991623535110900145190821631169159854457119691474665531075565e245"),
           e_float("-2.7797775581944616888811339285995580940563647902689752532435179719496567032130960448410109703312512547208827618533381963393101068707900966666212266516844634742195662630805236393333171982253827182353904982490756712308442249728563963037676304361696905652365588412371820899677568557790385553017846627202381008663456285299836675224300597943638626284894250351452563806575309697807623191271861379612069173786e249"),
           e_float("-6.8788259883240663399156650342089016170558849026768020985063430729950701807166263048806549131100487681046164144244327164699682484308127742252135520581879814803485038526539893537854827794126477511535435028253593896094884576416120857131953142531145432818659402109346097756723827649545622269063411631344096055994205219281254097334088338292329816837442940259729037384884544367674560490418264106457752939121e250"),
           e_float("-4.5892044305246331482578158054641449654728856869269287778583215449565999873301857648524093421831651945183722135763612425620592126750835098142590660824928735635411721871133875862638990759688521996313571556020480828388466041573352169031679110489049928368866984577992596764094822358206740969367035073728334036667363060321950224729089352524180553447791960409798103070540329003939714905063517843979018642141e251"),
           e_float("-1.2025281778894063599767946952202103439130239863982084323867191155689660366526308445902433353261585903442652854829938258331656846317560259521422504477694518903577099624327128571860451047163211781896703674969085380695691863056149164621456898226803432080767435251892400794836893081746772760782574285343412156977919251305179094812373516062359189075053993250658719856640503693808841889893485611369425330363e255"),
           e_float("-2.2835580731013169552708097940923373298584157887210805678643788505264407848762013389700773212616782534104705430325859216696318943940138168768284187232273753760902053437353494957544456125397868278859854626544956469827433013348223616727131984160233610612377913084588552463665653435984378948464138589866087153371199437735433438898549371369557812685971409736142473085560995208498611969626301885009436759931e262"),
           e_float("-6.4429776561216944182275665511791979728622238482738511530422565730305608723918958712375977748482456412359495437943993930513390393891795196120067226267916773702416710197896098559013963044549085074782394389675069347391570396962220049455597510969635738364632292489556808473436335673784144027595670747817585799254786108168110823905460464915597275178038541090834795072232628181099071302858431119213906582403e277"),
           e_float("-5.3968059437880761302883465403328870790668568847715760531522371408347198347691049521612717936074339856258248966001445247490361424664776589917450685607297626168493562174736853711545990850330989850451282721469860562606632356861373235010536577289031657971737073206673390001234078742476395869890030991037146200246556297505610159201868442995344524708232105598796123516489804493193307712821878236573585707427e373"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00685_legendre_pnu_vary_05(const bool b_write_output)
    {
      return TestCase_case_00685_legendre_pnu_vary_05().execute(b_write_output);
    }
  }
}