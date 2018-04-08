package test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

public class SplicingGraph {
	Map<Long, Long> used_kmers = new HashMap<Long, Long>();
	public class Node {
		private Vector<Long> parents = new Vector<Long>();
		private Vector<Long> children = new Vector<Long>();
		private String sequence;

		public boolean addParents(long parent) {
			if (parent < 0) {
				return false;
			} else {
				if (parents.size() != 0) {
					for (int i = 0; i < parents.size(); ++i) {
						if (parents.get(i) == parent) // if exist already
							return false;
					}
				}
				this.parents.add(parent);
				return true;
			}
		}

		public List getParents() {
			return parents;
		}

		public boolean addChildren(long child) {
			if (child < 0) {
				return false;
			} else {
				if (children.size() != 0) {
					for (int i = 0; i < children.size(); ++i) {
						if (children.get(i) == child) // if exist already
							return false;
					}
				}
				this.children.add(child);
				return true;
			}
		}

		public Vector getChildren() {
			return children;
		}

		public void setSequence(String seq) {
			sequence = seq;
		}

		public String getSequence() {
			return sequence;
		}

		public boolean isChild(long child) {
			if (child < 0) {
				return false;
			} else {
				for (int i = 0; i < children.size(); i++) {
					if (children.get(i) == child) {
						return true;
					}
				}
				return false;
			}
		}

		public boolean isParent(long parent) {
			if (parent < 0) {
				return false;
			} else {
				for (int i = 0; i < parents.size(); i++) {
					if (parents.get(i) == parent) {
						return true;
					}
				}
				return false;
			}
		}

		public boolean deleteChild(long child) {
			// TODO Auto-generated method stub
			if (child < 0) {
				return false;
			} else {
				for (int i = 0; i < children.size(); i++) {
					if (children.get(i) == child) {
						children.remove(i);
						return true;
					}
				}
				return false;
			}
		}

		public boolean deleteParent(long parent) {
			// TODO Auto-generated method stub
			if (parent < 0) {
				return false;
			} else {
				for (int i = 0; i < parents.size(); i++) {
					if (parents.get(i) == parent) {
						parents.remove(i);
						return true;
					}
				}
				return false;
			}
		}

		public void clearChildren() {
			// TODO Auto-generated method stub
			children.clear();
		}

		public void clearParents() {
			parents.clear();
		}

		public void clearAll() {
			// TODO Auto-generated method stub
			children.clear();
			parents.clear();
			sequence = "";
		}
	}

	public void add_node() {
		Node node = new Node();
		node.setSequence("zhangjianxi");
		System.out.println(node.getSequence());
	}
	
	private  String forward_extend(long kmer_int,Vector<Long> bifurcation,kmerHash kh){
		List<Map.Entry<Long,Long>> candidates=new ArrayList<Map.Entry<Long,Long>>();
		boolean flag = false;
		long sum=0l;
		String kmer_string=baseOptions.intvalToKmer(kmer_int, kh.kmer_length);
		while(true){
			candidates=kh.get_forward_candidates(kmer_int,kh.kmer_hash);
			//System.out.println("处理kmer:"+baseOptions.intvalToKmer(kmer_int, kh.kmer_length));
			if(candidates.size()==0){
				break;
			}
			int count=0;
			long candidate = 0l;
			long cov=0l;
			//System.out.println("其candidates为：");
			for (Map.Entry<Long, Long> mapping : candidates) {
				candidate=mapping.getKey();
				//System.out.println("candidate:"+baseOptions.intvalToKmer(candidate, kh.kmer_length)+"   read_count:"+mapping.getValue());
				if(!has_been_used(candidate)){//如果有未被使用的candidate
				//	System.out.println("选中："+baseOptions.intvalToKmer(candidate, kh.kmer_length));
					flag=true;
					cov=mapping.getValue();
					//如果还有其他未被使用的candidate
					if(count<candidates.size()-1){
						bifurcation.add(kmer_int);
					}
					break;
				}
				count++;
			//	System.out.println();
			}
			if(flag==true){
				used_kmers.put(candidate, cov);
				sum+=cov;
				int base_last_int=(int) (candidate&3l);
				char base_last_char=baseOptions.intToBase(base_last_int);
				kmer_string+=base_last_char;
				kmer_int=candidate;
			}
		}
		return kmer_string;
	}

	//反向扩展
	private  String reverse_extend(long kmer_int,Vector<Long> bifurcation,kmerHash kh){
		List<Map.Entry<Long,Long>> candidates=new ArrayList<Map.Entry<Long,Long>>();
		boolean flag = false;
		long sum=0l;
		String kmer_string=baseOptions.intvalToKmer(kmer_int, kh.kmer_length);
		while(true){
			candidates=kh.get_reverse_candidates(kmer_int,kh.kmer_hash);
			System.out.println("处理kmer:"+baseOptions.intvalToKmer(kmer_int, kh.kmer_length));
			if(candidates.size()==0){
				break;
			}
			int count=0;
			long candidate = 0l;
			long cov=0l;
			System.out.println("其candidates为：");
			for (Map.Entry<Long, Long> mapping : candidates) {
				candidate=mapping.getKey();
				System.out.println("candidate:"+baseOptions.intvalToKmer(candidate, kh.kmer_length)+"   read_count:"+mapping.getValue());
				if(!has_been_used(candidate)){//如果有未被使用的candidate
					System.out.println("选中："+baseOptions.intvalToKmer(candidate, kh.kmer_length));
					flag=true;
					cov=mapping.getValue();
					//如果还有其他未被使用的candidate
					if(count<candidates.size()-1){
						bifurcation.add(kmer_int);
					}
					break;
				}
				count++;
				System.out.println();
			}
			if(flag==true){
				used_kmers.put(candidate, cov);
				sum+=cov;
				int base_first_int=(int) (candidate>>(2*kh.kmer_length-2));
				char base_first_char=baseOptions.intToBase(base_first_int);
				kmer_string=base_first_char+kmer_string;
				kmer_int=candidate;
			}
		}
		return kmer_string;
	}

	
	private  boolean has_been_used(Long candidate) {
		if (used_kmers.containsKey(candidate)) {
			return true;
		} else {
			return false;
		}
	}

	public void init_graph(kmerHash kh) {
		// TODO Auto-generated method stub
//		System.out.println(kh.list.get(0));
//		System.out.println("开始kmer:"+baseOptions.intvalToKmer(kh.list.get(0).getKey(), kh.kmer_length));
//		System.out.println("**************");
//		used_kmers.put(kh.list.get(0).getKey(), kh.list.get(0).getValue());
//		System.out.println(used_kmers.toString());
//		System.out.println("结束！");
		Vector<Long> bifurcation=new Vector<Long>();
		String right=forward_extend(kh.list.get(0).getKey(),bifurcation,kh);
		String left=reverse_extend(kh.list.get(0).getKey(),bifurcation,kh);
		String trunk=left+right.substring(kh.kmer_length);
	//	System.out.println("左边扩展："+left);
	//	System.out.println("右边扩展："+right);
		System.out.println("******************"+trunk);
		System.out.println();
	}
}
