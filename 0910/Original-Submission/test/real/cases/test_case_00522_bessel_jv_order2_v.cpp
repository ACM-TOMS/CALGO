
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00522_bessel_jv_order2_v : public TestCaseReal
    {
    public:
      TestCase_case_00522_bessel_jv_order2_v() { }
      virtual ~TestCase_case_00522_bessel_jv_order2_v() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00522_bessel_jv_order2_v");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(11u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          const INT32 k2 = k * k;
          const INT32 k3 = k * k2;
          const INT32 k7 = (k2 * k2) * k3;
          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(ef::euler_gamma() + 100, k7 + ef::pi());
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("3.8746788152542443992435960739507927098627716421898715914816146904119845610204754595617173496541700433390939784084703153103375697332444155494021051637322294196262984679075474430843056297109105059364143625029953165309204759327006262258572497603817644819833693685965316729183291021226954369624542692260630358115311673531680734506279978798957348567330876364619278712856565502461227442128369611898342438709e-140"),
           e_float("4.4820348296116504424731311488540622073070183747619726772639305874485941832877751485120073995416534734790055181902155392526629344392660814105979647200801469903618182810672207671298905850392146159155118588596247122143680959564051556419896778510978005573920663057389425366728043624544279551854539972708551286761951599089474576081348169204283194239128680475003751677531661642919954231053354535971406434863e-128"),
           e_float("0.065184629679042834017307159859779273402953705114597634757217851861524126043344531286835520641550228773321392282479628492882562012427846374759294411420100578780955026979004730632374853320469403855766704441413703921942768304183477135922114701122111974405934169026692339827629254181602052976572487488683835714591345186483038012544228003375432977671703573678423338850474210906905519622421564005217356398482"),
           e_float("-0.0082025592342592231916165199222967263943977639505629239208920375038086588648915753062191494596016707613806372927477659911940261056312703355235035264248614885490607907766998848081955905747905585864355522509820138416813526007492481373948413990770933990155433770281466229919384301262044216552271216459284110819213801466745102183853702644213011101383791643279650845938286483088624571837652168402878427115512"),
           e_float("0.0043908642911451189056728054386486627607997718784923640703731450442947037774781897072723553556737479004930335221980174404880751697490374068791156055739923183332884776696583824729840393273642132852911360340792021687686337064500981879959005825456193983739450317350509429613181668858839035807041283658025865786076340142664411383032850016182191872609222673709726584331723470045048186245244976300709902392806"),
           e_float("0.000518523896175565808062066318719844925036421031205757103265527805570209801535642136975665609193597916520194798911779427651890104518098690416196938865009036471922256024477082555790737154842902489240305470442605529928894753927793153479810393860117433713842426586797942036365897519186000560506938081808052294664371473298676728091785453415497839751926388614340707017932881805198361566232277367920126963443"),
           e_float("-0.0013713640028541776838462224646658411690183612472166385075450063419652495955784235953529778873388553339399323926794154146288516704169177300191650240155851518439175299116665042837702326615170012129684812079284017358850371240416641469927821414452467975090884363292181314707580589722595415323060063528491538260280044861596020393831892503436993812233586551647226946041361204458225279865666181907974141149692"),
           e_float("0.00041885151524150386203941366248701909042803790203181761494740707115384056657033402904541455019545074636306342426081424845993644420812207085897338044161035134512442531026138466813161609001627619720014330480202992011875874849551912884956812497310376400024091609364721861348576179536497171498901625897971908684645196202644698341421907582090434534137739128140983578624907058876247350384666624963762186894967"),
           e_float("-0.00029022185366351816799399061059566947381442077917215366325908507926184177152172186567969260896898611106141495180065112863838424827776505060026547324705069137480576104829495871444828598474052762366583151493931486131765653361883505239490559024637083801397299765853170491912353202916296765507974313736079553860592880650823669062569003327918993152838360587922056571484622580363233455454731286529622190041509"),
           e_float("-0.00028091621963214663443782711794384515399832502474360249946589760052463904473954494178328501115161994219580210741515733127322524452651103942901546627127937883239980614011570900477094842195986323678707458159694012564815129171193324783167230958422443754173544259687284506624556694063074491590601829256565901585861047092530931303280217172385774255941332029149266907296263894210649807359942705084746603216503"),
           e_float("-0.00013291901286661523697217036436982746335287653582373190236856032528836357554768246300955641004955999681192992319444732207873308903188545136819522961386440330190503803637923196982073513587818279937886644787891202534101281236650455504261343628700850986167401469590035018903923702920263857541187311055853206375334256058021140367504103705812339222611591192632336993621452368426936487735769035330644405112207"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00522_bessel_jv_order2_v(const bool b_write_output)
    {
      return TestCase_case_00522_bessel_jv_order2_v().execute(b_write_output);
    }
  }
}