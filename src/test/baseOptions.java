package test;
import java.util.*;
import java.math.*;

public class baseOptions {
	long kmer_int_val; 
	static final int g_kmer_length=25;
	static char int_to_base[] = { 'A','C','G','T' };
	//ASCII码：A65 C67 G71 T84 a97 c99 g103 t116
	int base_to_int[] = { 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //   0-19
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  20-39
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, //  40-59
	255, 255, 255, 255, 255,   0, 255,   1, 255, 255, 255,  2 , 255, 255, 255, 255, 255, 255, 255, 255, //  60-79
	255, 255, 255, 255,   3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   0, 255,   1, //  80-99
	255, 255, 255,   2, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,   3, 255, 255, 255, // 100-119
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 120-139
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 140-159
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 160-179
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 180-209
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 200-219
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, // 220-239
	255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255                      // 240-255
	};

	boolean containNotAGCT(String kmer) {
		for (int i = 0;i < kmer.length();i++) {
			if (base_to_int[kmer.charAt(i)] > 3) {
				return true;
		    }
		}
		return false;
	}

	String revcomp(String kmer) {
		String reverseString="";
		for (int i = kmer.length() - 1;i >= 0;i--) {
			char tempc = kmer.charAt(i);
			switch (tempc)
			{
			case 'A':
				reverseString += 'T';
				break;
			case 'C':
				reverseString += 'G';
				break;
			case 'G':
				reverseString += 'C';
				break;
			case 'T':
				reverseString += 'A';
				break;
			case 'a':
				reverseString += 't';
				break;
			case 'c':
				reverseString += 'g';
				break;
			case 'g':
				reverseString += 'c';
				break;
			case 't':
				reverseString += 'a';
				break;
			default:
				break;
			}
		}
		return reverseString;
	}

	static char intToBase(int c) {
		if (c < 0 || c>3) {
			System.out.println("[警告！] int_val不在0-3内！");
			//exit(1);
			return 'Z';
		}
		//cout << int_to_base[int_val] << endl;
		return int_to_base[c];
	}

	static int baseToInt(char base_val) {
		switch (base_val) {
		case 'G':
		case 'g':
			return(2);
		case 'A':
		case 'a':
			return(0);
		case 'T':
		case 't':
			return(3);
		case 'C':
		case 'c':
			return(1);
		default:
			return(-1);
		}
	}

	

	public static String intvalToKmer(long int_val,int kmer_length) {
		String kmer = "";
		int c;
		char base;
		for (int i = kmer_length-1;i >=0;i--) {
			c = (int) (int_val & 3);
			base = intToBase(c);
			kmer = base + kmer;
			int_val = int_val >> 2;
		}
		return kmer;
	}

	public static float computeEntropy(String kmer) {
		int length = kmer.length();
		int a[] =new int[4];
		float entropy = 0;
		for (int i = 0;i < length;i++) {
			if(kmer.charAt(i)=='A'){
				a[0]++;
			}
			else if(kmer.charAt(i)=='C'){
				a[1]++;
			}
			else if(kmer.charAt(i)=='G'){
				a[2]++;
			}
			else if(kmer.charAt(i)=='T'){
				a[3]++;
			}
		}
		char base[] = { 'A','C','G','T' };
		int countl;
		double result;
		for (int i = 0;i < 4;i++) {
			countl=a[i];
			result = (float)countl / kmer.length();
			if (result > 0) {
				entropy += result * Math.log(1.0 / result) / Math.log(2.f);
			}
		}
		return entropy;
	}

	public static long kmerToIntval(String kmer) {
		long  intval = 0;
		for (int i = 0;i < kmer.length();i++) {
			int c=baseToInt(kmer.charAt(i));
			intval =intval<< 2;
			intval |= c;
			//System.out.println(intval);
		}
		
		return intval;
		//return 0;
	}

}
