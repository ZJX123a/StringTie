package test;
public class test {
	public static void main(String args[]){
		//System.out.println("hello");
		Node node=new Node();
		node.setSequence("vghgbkhkbhkbhj");
		System.out.println(node.getSequence());
		node.addParents(6l);
		System.out.println(node.addParents(6l));
		node.addParents(8l);
		node.addChildren(7l);
		node.addChildren(890l);
		System.out.println(node.getChildren());
		node.clearAll();
		System.out.println(node.getChildren());
		System.out.println(node.getParents());
		System.out.println(node.getSequence());
		
	}
}
