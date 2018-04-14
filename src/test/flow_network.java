package test;

public class flow_network {
public void init_graph(SplicingGraph sg,kmerHash kh) {
	// TODO Auto-generated method stub
	System.out.println(sg.node_set.size());
	//计算每个顶点的cov
	for(int i=0;i<sg.node_set.size();i++){
		int length=sg.node_set.get(i).getSequence().length();
		for(int j=0;j<length;j++){
			
		}
	}
}
}
