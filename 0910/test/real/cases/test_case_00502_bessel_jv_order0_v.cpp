
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00502_bessel_jv_order0_v : public TestCaseReal
    {
    public:
      TestCase_case_00502_bessel_jv_order0_v() { }
      virtual ~TestCase_case_00502_bessel_jv_order0_v() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00502_bessel_jv_order0_v");
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
          data[static_cast<std::size_t>(k)] = ef::cyl_bessel_j(ef::euler_gamma() + 1, k7 + ef::pi());
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 11u> a =
        {{
           e_float("0.46347831829332275532835785483293315761728103449813710193727418369438565630756681593491307576069365847277533545652295846590202983145728915618274900094625468796863591287242273079351662281716288365696999893358902291993624410115965549559728332709472254364615943413139465760515244096949301062504285970754061407147758147790668907458601687894546812550123830329756885725835133706069610479857308006264598342228"),
           e_float("0.1675569817734135365527010411994958447611596782839743495739097125185155078183750354123364129354539506782860017595332303359364772637771205032634163596710181819854870705595046466198117274004991128792551048604783324743437596223486464236074258660825070892949009613322255115157851225552693458251435264881514287019295009887737434177275708159134229991328036462806537207024307116008611665734629459263145189088"),
           e_float("-0.042318785252148125805210038409290909588836031247528875686047758864725806566665610375332811820077885137813421831519416659773696028147694488671272925244276653678333005445466894076791731558419443652869333654800926490310755407123389246430061620409053328897822012211621176778537882724764727176586693923947540967239599231609365223444220370927795531551222080024145985264487021432400051351160095470279984638078"),
           e_float("0.016125206850889549058500268502092075398910930766834915418135454244589597502109950853694849901840078544825478135014072341762664839753271096786436545616393825815475092514260619806391258724006190136162274346305973372655656092431251165159276107768980914233078342758124061790295467749343498886350806515541642698074423368270995147566769781437815702406023134882691687783504447468961461690913464589475118168535"),
           e_float("-0.0055482861103279078761886179272922276531399585283747747174903575369080297925825241216718122436548400866676529095035788505266882764235683109826727085132620637264072638944313382304564755194362065259059744406754012212119447058540210655426048920620161852042809512419949508772577049033378652162963894434773202239993533985489902786896234527623864800741868046532232268133799756671470819504246224243422464212793"),
           e_float("0.0027676375757478753484493432797287373857554309177736970368027488724523767139351386843294649349957443172766060331525446362984776738408209952622009999520643919109387604535681416228586855556009016172955103939604158278961256745714865015846051369942786688253062051185492426413589528537385163028748052689572093014958220756027384363854713005329190230371321102129020907336050651399447279152117561274413631957216"),
           e_float("0.00065196134338350308131553875398686872647874683231268014073309522496445110725391120395571331554666883171622764703078703472478922764744447751780560243482862201469147515612897068948059189769301089943674511370085949571828255762679078799459350037075918229608530704698092568647358967328639651074776428981564290154687011135718304111045572766583784460221307038820034004981244289747837629205664105007999469736582"),
           e_float("0.00077045052870018687843153729518179195638108185397446291830531504035150529367404787180763984547810224869010942352220955555911896042749441499054170324310597244914130549971620741786809455981005050369570713339577876494837502423254295379366609610368059174594550103264337445836143799665935215123276714432413886654350229330287640966302475104880214519157316034068012344699147052540908787269973202757317960261474"),
           e_float("0.00046903033523137889208482793253142201752540546030476501401693491151813735685177125352348565922168105592788981630032667775206728205792602978701553034124324626398078855066186249625454502815301212965206389314267027264675404183602780796534632584911306474552664003192617739708319902821911989393395851146160426566120387568543228850053154716671457906930237878078280403031230883957573842964124639207417279798894"),
           e_float("0.00023307866664188961980725088113687963083937818079977823458101518576997513530511447161396237953429460118031758626640123389298781561386081451711073229722617722490007829974193134566804589943917670384924627211467644520445374991752830967495765431451551722256397314357584050236969441080862631928598004425347288526749451687038235930120947361772485640034064732281020376879081270623000423420166108738130328869784"),
           e_float("-0.00021439602803610988753662811845981055899993121165044244170021537506922943144195178129870110322779457982313922470048887098507087733786399436263931278279932646978153535840694353523794722211817692931628850394142825341466964209099756267799396598750058099505290462492955706982312968986678164380388852087186453719992687523042917000193026990276788221453259006834568386801473865228643125233475614944509989790155"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00502_bessel_jv_order0_v(const bool b_write_output)
    {
      return TestCase_case_00502_bessel_jv_order0_v().execute(b_write_output);
    }
  }
}