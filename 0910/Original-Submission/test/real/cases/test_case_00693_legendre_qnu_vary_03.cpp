
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00693_legendre_qnu_vary_03 : public TestCaseReal
    {
    public:
      TestCase_case_00693_legendre_qnu_vary_03() { }
      virtual ~TestCase_case_00693_legendre_qnu_vary_03() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00693_legendre_qnu_vary_03");
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
          const INT32   n  = +((10 - k) * 13);
          const e_float u  = -((10 - k) * 11) - sqrt_1_5;
          data[k] = ef::legendre_q(n, u, x);
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("-3.1108017571383742354702104117732546946728782986060144903462008410103706020410804792093323655791551391126293225795906871982484424600894016095668697101233093768745116510172839403444377232998279180978502665959832622292881555190214851002546532808809102189988057010038565808013426992914735496069614001735059846029013853545333456059597548974716017752436767003983372674588859739447923714528461264163848487665e-136"),
           e_float("-3.4696809704424503071538234159891564239381068596900290993461498106098309977850360325107007243309459566902485906352022033244546315305338764519384786773480560698992296257719431373726723894050746815062027302306687810000862930700994658458919927688244611736978406075606944924025217138886921594985764714180555443124414209440414692085114694512462117807833081222552415552881704800691524629104898538052932371719e-191"),
           e_float("-7.9663275431354582641607154322894872811079245088842867842889690972819394000557706577858786820679390773749610469955118396740166627578042465217894550284326857589793391464427460656008824346831711040888753082221504068766889020610935224218150543805318271693479133594294038201114095387065218980140543920449599203101556990909420442993460401041393123361750843173843084441476268051867198085884771978704099897973e-174"),
           e_float("-6.978264232989513700811887178246164595179886496481729661518585532103881147561051453963744899549124782525843871495628724274879485421626063240780173306271620856261274251500117730669781556234531623805858629108073496159993442714348346759965366844204301497321462314578044385737109948593202962990199104466744621753035570612713345079820917986186863671702029874570747703269935124233326482476948998071692557738e-148"),
           e_float("-9.6199352909401953503700298855674143218753390588359620924779017598651176306692834480018147850315405354667497786522362313895596741378694746653788735303871649856342107958998361157693706602227371392728653518519233790979358555649739556452822105417385115156173754479035791923859993977172386125963106309145601149852982191760694117470195108129083540721480085942857722278704700966637650981834413401816569394296e-123"),
           e_float("-2.4968192495647831494382388995262650913737965044737372463063077517247798447971223540450890587083135162294730794949924829816279223035035948693776929739125690028552575836495898477656380033010998562633435873871274292096058647645092215103665036045062255222621531248147074067293771781587676362000787158489833116747997104351085894508688000247596201647804490900588660151064159215478623674140518726193208190805e-98"),
           e_float("-2.2455374457874457153660550177661570672010612410771921966959371015826126564755263682268675000781647054749265803912947049385015197482214409209737804607281177081962121726608572559051511481740083439835617695560274182404519271168590573678028918491996613394622177760661119583600327258829656221198703517646898550249609486197955584291037338163939269780031958010159077309053441589627062051556850885752092742654e-75"),
           e_float("-2.1790549247281909715957610243491402348443732692139280651781903193044378722538186489357532065400313324307123612346950602100289407106074807215231496589668497492473264962199304180330120065890192366222924960536169883633128326516105091126881581712845970520876494983671403966580212696997806953978792846342836568371852114037048590250516739095396739131775202844222381506906342833689485706888527666105019108889e-52"),
           e_float("3.2921074828439840397626470992719002641787534583644601677171966601706956832014576725902915829916632024333795686346543858967881297028540152415100217185922366486211860775256683125076469713807480904658767437645490146380667289382267564904823801968513972632747811760596746809911697052685715270713213506540206053841290705365765051757921556020209871194810585198625748851373940560204962980031156403929513596657e-31"),
           e_float("4.8777715570338539563213396926005511387446957702963525530292187819158373211364677752090327114109196277006230951070856691115989695394365842505161370588419947261452516899553667063951458644325297635159375524010970121639972656418972247794192128703952912478941285090209097166216267478452212070474141613158697358176992737774922156967679459493248383503635910880349649406577136622626663778630175084457311881798e-12"),
           e_float("6.5639431660463211515685953867679157293312877173900879378048961496976585430044156154496290514807594732096304979983453176576189333871875892166475221804085642937449759050096006135468903801283303280552247585621524690438426221671675986005454996699738049323516169408150830872562442669093603511317283375371349655981002946284708880958805720230347178113203082942633980056297645748977171180560060566277928663811"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00693_legendre_qnu_vary_03(const bool b_write_output)
    {
      return TestCase_case_00693_legendre_qnu_vary_03().execute(b_write_output);
    }
  }
}