package bioproject;

import static bioproject.Five.ANS;
import static bioproject.Four.codon;
import static bioproject.Four.engname;
import static bioproject.Four.shname;
import static bioproject.Four.thname;
import static bioproject.Four.yor;
import java.util.*;

/**
 *
 * @author Nize
 */
public class Algorithms {

    private int c;

    protected String Reverse(String a) {
        return PrintFormat(new StringBuffer(a).reverse().toString());
    }

    protected String Complementary(String a) {
        String ans = "";
        for (int i = 0; i < a.length(); i++) {
            if ((a.charAt(i) + "").equalsIgnoreCase("A")) {
                ans += "T";
            } else if ((a.charAt(i) + "").equalsIgnoreCase("T") || (a.charAt(i) + "").equalsIgnoreCase("U")) {
                ans += "A";
            } else if ((a.charAt(i) + "").equalsIgnoreCase("G")) {
                ans += "C";
            } else if ((a.charAt(i) + "").equalsIgnoreCase("C")) {
                ans += "G";
            }
        }
        return PrintFormat(ans);
    }

    protected boolean CheckError(String a) {
        for (int i = 0; i < a.length(); i++) {
            if ((a.charAt(i) + "").equalsIgnoreCase("A")) {
            } else if ((a.charAt(i) + "").equalsIgnoreCase("T")) {
            } else if ((a.charAt(i) + "").equalsIgnoreCase("G")) {
            } else if ((a.charAt(i) + "").equalsIgnoreCase("C")) {
            } else if ((a.charAt(i) + "").equalsIgnoreCase("U")) {
            } else {
                return true;
            }
        }
        return false;
    }

    protected String DNATomRNA(String a) {
        char tmp[] = a.toCharArray();
        for (int i = 0; i < tmp.length; i++) {
            if ((tmp[i] + "").equalsIgnoreCase("T")) {
                tmp[i] = 'U';
            }
        }
        return PrintFormat(Arrays.toString(tmp).replaceAll("\\[|\\]|,|\\s", ""));
    }

    protected String mRNAToDNA(String a) {
        char tmp[] = a.toCharArray();
        for (int i = 0; i < tmp.length; i++) {
            if ((tmp[i] + "").equalsIgnoreCase("U")) {
                tmp[i] = 'T';
            }
        }
        return PrintFormat(Arrays.toString(tmp).replaceAll("\\[|\\]|,|\\s", ""));
    }

    protected char EqualCodon(String a) {
        switch (a) {
            case "UUU":
            case "UUC":
                //1
                return 'F';
            case "UUA":
            case "UUG":
                return 'L';
            case "UCU":
            case "UCC":
            case "UCA":
            case "UCG":
                return 'S';
            case "UAU":
            case "UAC":
                return 'Y';
            case "UAA":
            case "UAG":
            case "UGA":
                return '#';
            case "UGU":
            case "UGC":
                return 'C';
            case "UGG":
                return 'W';
            case "CUU":
            case "CUC":
            case "CUA":
            case "CUG":
                //2
                return 'L';
            case "CCU":
            case "CCC":
            case "CCA":
            case "CCG":
                return 'P';
            case "CAU":
            case "CAC":
                return 'H';
            case "CAA":
            case "CAG":
                return 'Q';
            case "CGU":
            case "CGC":
            case "CGA":
            case "CGG":
                return 'R';
            case "AUU":
            case "AUC":
            case "AUA":
                //3
                return 'I';
            case "AUG":
                return '*';
            case "ACU":
            case "ACC":
            case "ACA":
            case "ACG":
                return 'T';
            case "AAU":
            case "AAC":
                return 'N';
            case "AAA":
            case "AAG":
                return 'K';
            case "AGU":
            case "AGC":
                return 'S';
            case "AGA":
            case "AGG":
                return 'R';
            case "GUU":
            case "GUC":
            case "GUA":
            case "GUG":
                //4
                return 'V';
            case "GCU":
            case "GCC":
            case "GCA":
            case "GCG":
                return 'A';
            case "GAU":
            case "GAC":
                return 'D';
            case "GAA":
            case "GAG":
                return 'E';
            case "GGU":
            case "GGC":
            case "GGA":
            case "GGG":
                return 'G';
            default:
                break;
        }
        return '!';
    }

    protected String ConpareAndPlus(String tmpcodon) {
        String ans = "";
        char eqa = EqualCodon(tmpcodon);
        switch (eqa) {
            case '!':
                ans += "";
                break;
            case '*':
                ans += "Met ";
                c++;
                break;
            case '#':
                ans += "Stop ";
                break;
            case 'F':
                ans += "F ";
                c++;
                break;
            case 'L':
                ans += "L ";
                c++;
                break;
            case 'S':
                ans += "S ";
                c++;
                break;
            case 'Y':
                ans += "Y ";
                c++;
                break;
            case 'C':
                ans += "C ";
                c++;
                break;
            case 'W':
                ans += "W ";
                c++;
                break;
            case 'P':
                ans += "P ";
                c++;
                break;
            case 'H':
                ans += "H ";
                c++;
                break;
            case 'Q':
                ans += "Q ";
                c++;
                break;
            case 'R':
                ans += "R ";
                c++;
                break;
            case 'I':
                ans += "I ";
                c++;
                break;
            case 'T':
                ans += "T ";
                c++;
                break;
            case 'N':
                ans += "N ";
                c++;
                break;
            case 'K':
                ans += "K ";
                c++;
                break;
            case 'V':
                ans += "V ";
                c++;
                break;
            case 'A':
                ans += "A ";
                c++;
                break;
            case 'D':
                ans += "D ";
                c++;
                break;
            case 'E':
                ans += "E ";
                c++;
                break;
            case 'G':
                ans += "G ";
                c++;
                break;
            default:
                break;
        }
        return ans;

    }

    protected String ChangeDNAToAminoAcid(String a) {
        int len = a.length();
        String tmpcodon, answer = "";
        c = 0;
        a = DNATomRNA(a).toUpperCase().replaceAll("\n", "");             //change DNA to mRNA     
        ANS += ("> Sequence convert 1 frame \n");
        for (int i = 0; i < len; i += 3) {
            if (i + 2 < len) {
                answer += ConpareAndPlus("" + a.charAt(i) + a.charAt(i + 1) + a.charAt(i + 2));
            }
        }
        String p1[] = answer.split(" ");
        for (int i = 1; i <= p1.length; i++) {
            if (i % 11 == 0) {
                ANS += "\n";
            }
            ANS += p1[i - 1] + "  ";
        }
        ANS += ("\nจำนวนกรดอะมิโนที่ได้จากการแปลงรหัส DNA = " + c + " ตัว\n");
        c = 0;
        answer = "";

        ANS += ("\n> Sequence convert 2 frame \n");
        for (int i = 1; i < len; i += 3) {
            if (i + 2 < len) {
                answer += ConpareAndPlus("" + a.charAt(i) + a.charAt(i + 1) + a.charAt(i + 2));
            }
        }
        String p2[] = answer.split(" ");
        for (int i = 1; i <= p2.length; i++) {
            if (i % 11 == 0) {
                ANS += "\n";
            }
            ANS += p2[i - 1] + "  ";
        }
        ANS += ("\nจำนวนกรดอะมิโนที่ได้จากการแปลงรหัส DNA = " + c + " ตัว\n");
        c = 0;
        answer = "";

        ANS += ("\n> Sequence convert 3 frame \n");
        for (int i = 2; i < len; i += 3) {
            if (i + 2 < len) {
                answer += ConpareAndPlus("" + a.charAt(i) + a.charAt(i + 1) + a.charAt(i + 2));
            }
        }
        String p3[] = answer.split(" ");
        for (int i = 1; i <= p3.length; i++) {
            if (i % 11 == 0) {
                ANS += "\n";
            }
            ANS += p3[i - 1] + "  ";
        }
        ANS += ("\nจำนวนกรดอะมิโนที่ได้จากการแปลงรหัส DNA = " + c + " ตัว\n");
        c = 0;
        answer = "";

        //set reverse and complement
        a = Complementary(new StringBuffer(a).reverse().toString()).replaceAll("\n", "");
        a = DNATomRNA(a).replaceAll("\n", "");

        ANS += ("\n> Sequence convert 4 frame \n");
        for (int i = 0; i < len; i += 3) {
            if (i + 2 < len) {
                answer += ConpareAndPlus("" + a.charAt(i) + a.charAt(i + 1) + a.charAt(i + 2));
            }

        }
        String p4[] = answer.split(" ");
        for (int i = 1; i <= p4.length; i++) {
            if (i % 11 == 0) {
                ANS += "\n";
            }
            ANS += p4[i - 1] + "  ";
        }
        ANS += ("\nจำนวนกรดอะมิโนที่ได้จากการแปลงรหัส DNA = " + c + " ตัว\n");
        c = 0;
        answer = "";

        ANS += ("\n> Sequence convert 5 frame \n");
        for (int i = 1; i < len; i += 3) {
            if (i + 2 < len) {
                answer += ConpareAndPlus("" + a.charAt(i) + a.charAt(i + 1) + a.charAt(i + 2));
            }
        }
        String p5[] = answer.split(" ");
        for (int i = 1; i <= p5.length; i++) {
            if (i % 11 == 0) {
                ANS += "\n";
            }
            ANS += p5[i - 1] + "  ";
        }
        ANS += ("\nจำนวนกรดอะมิโนที่ได้จากการแปลงรหัส DNA = " + c + " ตัว\n");
        c = 0;
        answer = "";
        ANS += ("\n> Sequence convert 6 frame \n");
        for (int i = 2; i < len; i += 3) {
            if (i + 2 < len) {
                answer += ConpareAndPlus("" + a.charAt(i) + a.charAt(i + 1) + a.charAt(i + 2));
            }
        }
        String p6[] = answer.split(" ");
        for (int i = 1; i <= p6.length; i++) {
            if (i % 11 == 0) {
                ANS += "\n";
            }
            ANS += p6[i - 1] + "  ";
        }
        ANS += ("\nจำนวนกรดอะมิโนที่ได้จากการแปลงรหัส DNA = " + c + " ตัว\n");

        return ANS;
    }

    protected String GenePrediction(String a) {
        String tmpcodon, tmp1, an, re = "";
        ArrayList<String> t1 = new ArrayList<>();
        for (int i = 0; i < a.length(); i += 3) {
            if (i + 2 < a.length()) {
                tmpcodon = "" + a.charAt(i) + a.charAt(i + 1) + a.charAt(i + 2);
                if (tmpcodon.equals("AUG")) {
                    an = tmpcodon;
                    if (i + 3 < a.length()) {
                        for (int j = i + 3; j < a.length(); j += 3) {
                            if (j + 2 < a.length()) {
                                tmpcodon = "" + a.charAt(j) + a.charAt(j + 1) + a.charAt(j + 2);
                                if (tmpcodon.equals("UAA") || tmpcodon.equals("UAG") || tmpcodon.equals("UGA")) {
                                    an += tmpcodon;
                                    break;
                                } else {
                                    an += tmpcodon;
                                }
                            }
                        }
                    }
                    tmpcodon = "" + an.charAt(an.length() - 3) + an.charAt(an.length() - 2) + an.charAt(an.length() - 1);
                    if (tmpcodon.equals("UAA") || tmpcodon.equals("UAG") || tmpcodon.equals("UGA")) {
                        t1.add(mRNAToDNA(an));
                    }
                }
            }
        }
        for (int i = 0; i < t1.size(); i++) {
            re += "ยีน #" + (i + 1) + "\n";
            re += PrintFormat(t1.get(i).replaceAll("\\[|\\]|,|\\s", "")) + "\n\n";
        }
        return re;
    }

    protected String PrintFormat(String a) {
        String ans = "";
        for (int i = 1; i <= a.length(); i++) {
            ans += a.charAt(i - 1);
            if (i % 30 == 0) {
                ans += "\n";
            }
        }
        return ans;
    }

    protected String CutSmall(String a) {
        String ans = "";
        for (int i = 0; i < a.length(); i++) {
            if (a.charAt(i) >= 'A' && a.charAt(i) <= 'Z') {
                ans += a.charAt(i);
            }
        }
        return PrintFormat(ans);
    }

    protected String CutCapital(String a) {
        String ans = "";
        for (int i = 0; i < a.length(); i++) {
            if (a.charAt(i) >= 'a' && a.charAt(i) <= 'z') {
                ans += a.charAt(i);
            }
        }
        return PrintFormat(ans);
    }

    protected String ChangeBigToSmall(int a, int b, String str) {
        String ans = "", ans1 = "";
        if (a <= str.length() && b <= str.length() && a <= b) {
            for (int i = 0; i < a - 1; i++) {
                ans += str.charAt(i);
            }
            for (int i = a - 1; i < b; i++) {
                ans1 += str.charAt(i);
            }
            ans += ans1.toLowerCase();

            for (int i = b; i < str.length(); i++) {
                ans += str.charAt(i);
            }
        } else {
            ans = "*";
        }
        return PrintFormat(ans);
    }

    protected String ChangeSmallToBig(int a, int b, String str) {
        String ans = "", ans1 = "";
        if (a <= str.length() && b <= str.length() && a <= b) {
            for (int i = 0; i < a - 1; i++) {
                ans += str.charAt(i);
            }
            for (int i = a - 1; i < b; i++) {
                ans1 += str.charAt(i);
            }
            ans += ans1.toUpperCase();

            for (int i = b; i < str.length(); i++) {
                ans += str.charAt(i);
            }
        } else {
            ans = "*";
        }
        return PrintFormat(ans);
    }

    protected int Search(String a) {
        int ans = -1;
        for (int i = 0; i < 64; i++) {
            if (a.equalsIgnoreCase(codon[i]) || a.equalsIgnoreCase(shname[i])
                    || a.equalsIgnoreCase(yor[i] + "") || a.equalsIgnoreCase(engname[i])
                    || a.equalsIgnoreCase(thname[i])) {
                return i;
            }
        }
        return ans;
    }

    protected String SameCodon(int tmp) {
        String name = engname[tmp], ans = "";
        String tm[] = name.split(" ");
        name = tm[0];
        for (int i = 0; i < 64; i++) {
            tm = engname[i].split(" ");
            if (name.equalsIgnoreCase(tm[0])) {
                ans += " " + codon[i] + ",";
            }
        }
        if (ans.length() > 0) {
            ans = ans.substring(1, ans.length() - 1);
        }
        return ans;
    }
}
