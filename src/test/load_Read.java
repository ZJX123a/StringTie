package test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.Vector;

public class load_Read {
	public static Vector<String> read_vector;
	   public static  void load_reads() throws IOException{
		   File file=new File("sim50bp_2.fa");
		   BufferedReader reader=null;
		   read_vector=new Vector<String>();
		   try {
			reader=new BufferedReader(new FileReader(file));
			String tempString=null;
			int length=0;
			while((tempString=reader.readLine())!=null){
				tempString=reader.readLine();
				length++;
				read_vector.addElement(tempString);
				//System.out.println(length+"     "+tempString);
			}
			/*Iterator iterator=read_vector.iterator();
			while(iterator.hasNext()){
				System.out.println(iterator.next());
			}
			System.out.println(length);*/
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	   }
public static void main(String args[]) throws IOException{
	load_reads();
	kmerHash kh=new kmerHash();
	kh.readsToKmer();
}


}