package test;

import java.io.IOException;
import test.create_graph;
import java.util.*;
import test.baseOptions;
import test.load_Read;

public class kmerHash {
	final int kmer_length = 25;
	Map<Long, Vector<Integer>> kmer_hash = new HashMap<Long, Vector<Integer>>();
	float g_min_ratio_non_error = 0.5f;
	Map<Long, Long> kmer_map = new HashMap<Long, Long>();
	List<Map.Entry<Long, Long>> list;

	public void readsToKmer(Map<Long, Vector<Integer>> kmer_hash) {
		// TODO Auto-generated method stub
		String temp_kmer;
		long intval = 0;
		for (int i = 0; i < load_Read.read_vector.size(); i++) {
			String read = load_Read.read_vector.get(i);
			for (int j = 0; j < read.length() - kmer_length + 1; j++) {
				temp_kmer = read.substring(j, j + kmer_length);
				intval = baseOptions.kmerToIntval(temp_kmer);
				Vector<Integer> readset = new Vector<Integer>();
				if (kmer_hash.containsKey(intval)) {
					readset = get_readset(intval, kmer_hash);
					readset.add(i);
				} else {
					readset.add(i);
				}
				kmer_hash.put(intval, readset);
			}

		}
	}

	private Vector<Integer> get_readset(long intval, Map<Long, Vector<Integer>> kmer_hash) {
		return (Vector<Integer>) kmer_hash.get(intval);
	}

	private boolean ifexist(Map<Long, Vector<Integer>> kmer_hash, long intval) {
		if (kmer_hash.containsKey(intval)) {
			return true;
		} else {
			return false;
		}
	}

	private long get_readset_count(Map<Long, Vector<Integer>> kmer_hash, long intval) {
		if (ifexist(kmer_hash, intval)) {
			Vector readset = (Vector) kmer_hash.get(intval);
			long count = readset.size();
			return count;
		} else {
			return 0;
		}
	}

	private void delete_bad_kmers(Map<Long, Vector<Integer>> kmer_hash) {
		Vector delete_list = new Vector();
		for (Iterator<Long> iterator = kmer_hash.keySet().iterator(); iterator.hasNext();) {
			List<Map.Entry<Long,Long>> candidates=new ArrayList<Map.Entry<Long,Long>>();
			Long key = iterator.next();
			Long candidate;
			candidates = get_forward_candidates(key, kmer_hash);
			long dominate_count = 0;
			if (candidates.size() > 0) {
				for (Map.Entry<Long, Long> mapping : candidates) {
					candidate=mapping.getKey();
					if (mapping.getValue() > 0) {
						long candidate_count = mapping.getValue();
						if (dominate_count == 0) {
							dominate_count = candidate_count;
						} else if ((float) candidate_count / dominate_count < g_min_ratio_non_error) {
							delete_list.add(candidate);
						}
					}
				}
			}
		}
		 System.out.println("之前的hash");
		 System.out.println(kmer_hash.size());
		if (!delete_list.isEmpty()) {
			for (int j = 0; j < delete_list.size(); j++) {
				long delete_kmer = (long) delete_list.get(j);
				kmer_hash.remove(delete_kmer);
			}
		}
System.out.println(kmer_hash.size());
	}

	private boolean if_find_kmer(long intval, Map<Long, Vector<Integer>> kmer_hash) {
		if (kmer_hash.containsKey(intval)) {
			return true;
		} else {
			return false;
		}
	}

	public List get_forward_candidates(long seed_kmer, Map<Long, Vector<Integer>> kmer_hash) {
		if (ifexist(kmer_hash, seed_kmer)) {
			Map forward_candidates = new HashMap();
			List<Map.Entry<Long, Long>> list_forward = new ArrayList<Map.Entry<Long, Long>>();
			long temp_intval = seed_kmer << (65 - kmer_length * 2) >> (63 - kmer_length * 2);
			long a[] = { 0l, 1l, 2l, 3l };
			for (int i = 0; i < 4; i++) {
				long temp2 = temp_intval;
				temp2 |= a[i];
				temp2 = baseOptions.kmerToIntval(baseOptions.intvalToKmer(temp2, kmer_length));
				if (ifexist(kmer_hash, temp2)) {
					long read_count = get_readset_count(kmer_hash, temp2);
					forward_candidates.put(temp2, read_count);
				}
				list_forward = new ArrayList<Map.Entry<Long, Long>>(forward_candidates.entrySet());

				// 通过比较器实现比较排序
				Collections.sort(list_forward, new Comparator<Map.Entry<Long, Long>>() {
					public int compare(Map.Entry<Long, Long> o1, Map.Entry<Long, Long> o2) {
						return o2.getValue().compareTo(o1.getValue()); // 倒序
					}
				});
			}
			return list_forward;
		}
		return null;

	}

	public List get_reverse_candidates(long seed_kmer, Map<Long, Vector<Integer>> kmer_hash) {
		if (ifexist(kmer_hash, seed_kmer)) {
			Map reverse_candidates = new HashMap();
			List<Map.Entry<Long, Long>> list_reverse = new ArrayList<Map.Entry<Long, Long>>();
			long temp_intval = seed_kmer >> 2;
			long a[] = { 0l, 1l, 2l, 3l };
			for (int i = 0; i < 4; i++) {
				long temp2 = (a[i] << (2 * kmer_length - 2)) | temp_intval;
				temp2 = baseOptions.kmerToIntval(baseOptions.intvalToKmer(temp2, kmer_length));
				if (ifexist(kmer_hash, temp2)) {
					long read_count = get_readset_count(kmer_hash, temp2);
					reverse_candidates.put(temp2, read_count);
				}

				list_reverse = new ArrayList<Map.Entry<Long, Long>>(reverse_candidates.entrySet());

				// 通过比较器实现比较排序
				Collections.sort(list_reverse, new Comparator<Map.Entry<Long, Long>>() {
					public int compare(Map.Entry<Long, Long> o1, Map.Entry<Long, Long> o2) {
						return o2.getValue().compareTo(o1.getValue()); // 倒序
					}
				});
			}
			return list_reverse;
		} else {
			return null;
		}

	}

	private void print_kmerhash(Map<Long, Vector<Integer>> kmer_hash) {
		System.out.println(kmer_hash.size());
		for (Iterator<Long> iterator = kmer_hash.keySet().iterator(); iterator.hasNext();) {
			Long key = iterator.next();
			System.out.println("key-----" + baseOptions.intvalToKmer(key, kmer_length) + "     " + key);
			System.out.println("value--------" + kmer_hash.get(key));
		}
	}

	public List sort_kmer(Map<Long, Vector<Integer>> kmer_hash) {

		Long key;
		Long count;
		for (Iterator<Long> iterator = kmer_hash.keySet().iterator(); iterator.hasNext();) {
			key = iterator.next();
			count = (long) kmer_hash.get(key).size();
			kmer_map.put(key, count);
		}
		list = new ArrayList<Map.Entry<Long, Long>>(kmer_map.entrySet());

		// 通过比较器实现比较排序
		Collections.sort(list, new Comparator<Map.Entry<Long, Long>>() {
			public int compare(Map.Entry<Long, Long> o1, Map.Entry<Long, Long> o2) {
				return o2.getValue().compareTo(o1.getValue()); // 倒序
			}
		});

		// for (Map.Entry<Long, Long> mapping : list) {
		// System.out.println(baseOptions.intvalToKmer(mapping.getKey(),
		// kmer_length) + ":" + mapping.getValue()
		// + " intval:" + mapping.getKey());
		// }
return list;
	}

	public static void main(String args[]) throws IOException {
		load_Read.load_reads();
		kmerHash kh = new kmerHash();
		kh.readsToKmer(kh.kmer_hash);
		// kh.print_kmerhash(kh.kmer_hash);
		System.out.println("***********************");

		// print_kmerhash();
		kh.delete_bad_kmers(kh.kmer_hash);
		// kh.get_reverse_candidates(110653935488931, kh.kmer_hash);
		kh.sort_kmer(kh.kmer_hash);
		List listK=kh.sort_kmer(kh.kmer_hash);
		if(listK.size()==0){
			System.out.println("没有数据！");
			return;
		}
		System.out.println("***********************");
		// System.out.println(kh.list.toString());
		//create_graph cg = new create_graph();
		//cg.init_graph(kh);
		SplicingGraph sg=new SplicingGraph();
		sg.init_graph(kh);
		System.out.println("---------------------------------");
		// String kmer= baseOptions.intvalToKmer(110653935488931l,
		// kh.kmer_length);
		// long kmer_int=110653935488931l>>2;
		// kmer_int=kmer_int|(2l<<(2*kh.kmer_length-2));
		// System.out.println(baseOptions.intvalToKmer(kmer_int,
		// kh.kmer_length));
	}
}
