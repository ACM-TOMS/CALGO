
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00492_bessel_jv_large_x : public TestCaseReal
    {
    public:
      TestCase_case_00492_bessel_jv_large_x() { }
      virtual ~TestCase_case_00492_bessel_jv_large_x() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00492_bessel_jv_large_x");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(100 + ef::third(), ef::euler_gamma() + ef::million() * (k + 1));
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("0.00037302754764928617254382233977661587907175350354526083376286084788393869781758299287621371466372174250703574016545466925685858991483082663211467602253632389256014484699164469558563754127439960350140365392038423192084834782530040076957919647716870284901902287190119050841164846406557420358218712730181645859438287017420589571090931529525891063101764847158578066203209981003301869542885875058442721762945"),
           e_float("0.000071125683938506233717796205978423585925765218703061753107464858401057824660191989207770171116901113643597374084899590477360373278256102607320382921037227349126207325282707106620460962814228791703261675971840234170450182180918326811217442170367445246072582875494763008732087253811972816646272873669064612413310407130462034582009838976613378366351480844537236003834205834193474550548891694970145460745996"),
           e_float("-0.00010591656570048841191055492305152146570391339539627835506111424427136258112665420204980999700875778244255616641534623130980672771085347614044687970820735282207184076991882818832195771544209667179078205771647580894255714205784174189750383007840407873824357051149245133667877731581191870124204424362766808978416831181754773861868288239528512480232836549296300595237978552514104089296212047943510736126917"),
           e_float("-0.00022195038314182517805334289059667992905661584032735226204355754012071324646066660564296255721865852919407993492875951904225844520258957896081542875273686245817170779543578959592067450881079190684349452180492017566574315407121593704197551013575868980829628630627869301862819362245438611595014093486688408093421879101943766782754442303192840110944698031758037812324958181667305242447576493486769786774817"),
           e_float("-0.00028978931414086363594972180773344788979018593677878721627812037421362666498917489054747780040386034289774010949475736101078717724521536118246391339799957805061262813653912640895507320220830980602723926186952613366296685144620609641292687995864891934502426560648203075347470429511880877864432522958425285242541977179829925809733507138670013652231363405036234950962740014014010763249974628623076443467696"),
           e_float("-0.00031434175763802002960482396378563883886675400081062020299948296573386012964013917301751369841837394357569659885678457535969742987382207156990858651453809774152956920184385436282537380984272225623994708127074424875220375229569552966334787466355639872208553518151414057285193007401423269110043997610707842875329070875548795183563508407621151722084246212881639692989878650288837752350743494626285299538962"),
           e_float("-0.00030028502425553690356591359664042905223367864804618787348325416990330528275034140396606844755435890504015844362908377315745203640845870549922594697834378980215085888223852874010488726890250551176204856165815724558803792746991298485730812991863292813831360957454047670974137708876257283664547422917097829317932745247450078309669348167521794582542848147307100389350814544768205438194396933793312127498419"),
           e_float("-0.00025400249543710244155175855223172425983629175152366003895424914197056581799167677054594132654822373651598921468323535554484735845056333744383718320942564510776801568355589290259203527710554795437833456057814183835393133290905958299994560405057501629236416789506540981380383196020320909307248991083223635824132063099036857733286249962870733054125714784404086954516830460337130541017722619975355463929661"),
           e_float("-0.00018382129390094493823368691939558702701121438304065293731028519681219272898806369359709887657529879972127092426070135650623368018318316411895289052683438384146072615074190470332088605582965321872942211465248417766643753336395922379558151877950303867792713754142373939684321687938498487157990410702339723514107217194484453351584051733738707867977424882828811715014481463178416063090202503155568118154527"),
           e_float("-0.000099524986837897779133548063342984558212146897065281492323611664984349642568212000490439861690183404730799836570064278791138008900793705176266217252247561478061642741214198019134363036115099024470166899247697871238344505783127145072394311526500320154462416527964408107026284894767764379624042290137122074250446408670505953155192992817736250564934685223170943975202515057103562092251003686910171107882794"),
           e_float("-0.000011509209445936939210190220468780760159785520292870910635466879068158511552042663714131764000636609524980961817350060793940972092688019261902853939739944734007033125026848990143008756703741486041299491405464841840646397188142316686527162500686483985334771941353489747487912698497170685197046755551273784734955492432930951825296058511426301424945523375063209775631678628844208968200100058492380698060688"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00492_bessel_jv_large_x(const bool b_write_output)
    {
      return TestCase_case_00492_bessel_jv_large_x().execute(b_write_output);
    }
  }
}
