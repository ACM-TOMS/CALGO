
// Automatically generated file
#include <functions/functions.h>
#include <test/real/test_case_real.h>

namespace test
{
  namespace real
  {
    class TestCase_case_00957_polylog_x_large_neg_n : public TestCaseReal
    {
    public:
      TestCase_case_00957_polylog_x_large_neg_n() { }
      virtual ~TestCase_case_00957_polylog_x_large_neg_n() { }
    private:
      virtual const std::string& name(void) const
      {
        static const std::string str("TestCase_case_00957_polylog_x_large_neg_n");
        return str;
      }
      virtual void e_float_test(std::vector<e_float>& data) const
      {
        data.resize(6u);
        for(INT32 k = static_cast<INT32>(0); k < static_cast<INT32>(data.size()); k++)
        {
          data[static_cast<std::size_t>(k)] = ef::poly_logarithm(-500 -(k * 100), ef::five_k());
        }
      }
      virtual const std::vector<e_float>& control_data(void) const
      {
        static const std::tr1::array<e_float, 6u> a =
        {{
           e_float("-1.018269911957888776381477189131070889733776867764340225580056168077681396083445456488159085059570000097760306553300591678796549652273466716301732900549842593477699826875205424983381259712696954052811904764932131411015251939408019263952137069984069870541211582193290159535238956552153067143371474635367897453522714366218203443059413234236125441474586363559174561331156362910027113940167413589709063519e668"),
           e_float("-9.8648734279096679264197852611702344949210748619664445758109359853250967262692609257148439441619961766733949203003791262773604413064916274637795427535952736902399132759007783386195486178691379403587107957368930165383527349714135032982052571961565664555627676605812741230403544381881352659623501223789505208908362593078147384576932957009416234095690965418265176189692438211537190226003949049167375503127e848"),
           e_float("-1.7633400013683860058398452677554148568587286774054255714392271131108035174495996508119343337196198806374220110298219188083515159324890558315358601458103776045720010624526966505731883070264122212190510210101065456531012826957497268565211494202557489187194126905720449990430268753019965016294096062121215074683291755766553202124922745671422934500006811120690178091757249280761989070647849617523544035027e1037"),
           e_float("-5.2431120431767515890853112248269339953444212148044555737254486177527289612630149313372660218618767773291606003260314624776544092886097854518798067477089340810710392279098829861566850325657664679127249401943641761037590956705896770007915431414057080468748045837783807608108809719875716672806299248565837043469560150414565597146864237347370948978514939851714001868645137591218088356179365940299792134262e1231"),
           e_float("-4.2887552192155971484800749162663967875269920828875137146496886508184677916147654461117517717636181457017923537153852822387641207622990059996322902044898760074988771899500642444041568291097444393846704260601230136221541058351878254122119901477401519967667893206530718652610028330548247543737348625639217151448867895722528923065062655872774181003559970775197514036383435660630798563603177342367840703297e1431"),
           e_float("-2.3869838068477271810147082734420976196545562593378027304400914001261759737852568305996234047228758900914191410751589085195747143107497765619381848146004192552479756508849046228740278145347476302566101039968381409664166880410188852117177115742232415171037030764802409924345949035286850633346798113577190781131294910655696417142503609322140963776968886538897514260939746194017796307895046516764191150966e1636"),
        }};
        static const std::vector<e_float> v(a.begin(), a.end());
        return v;
      }
    };

    bool test_case_00957_polylog_x_large_neg_n(const bool b_write_output)
    {
      return TestCase_case_00957_polylog_x_large_neg_n().execute(b_write_output);
    }
  }
}
