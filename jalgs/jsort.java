/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class jsort{
	// here we have sorting routines. The ones based on the java collections
	// framework are more robust but perhaps a bit slower

	// this basic routine is adapted from Sedgewick 1978
	public void quicksort(float[] arr){
		int i,j,k;
		int ir=arr.length;
		int l=1;
		int jstack=0;
		int[] istack;
		float a,temp2;
		int M=7;

		istack=new int[50];
		for(;;){
			if(ir-l<M){
				for(j=l+1;j<=ir;j++){
					a=arr[j-1];
					for(i=j-1;i>=1;i--){
						if(arr[i-1]<=a)
							break;
						arr[i]=arr[i-1];
					}
					arr[i]=a;
				}
				if(jstack==0)
					break;
				ir=istack[(jstack--)-1];
				l=istack[(jstack--)-1];
			}else{
				k=(l+ir)>>1;
				temp2=0.0f;
				temp2=arr[k-1];
				arr[k-1]=arr[l];
				arr[l]=temp2;
				if(arr[l]>arr[ir-1]){
					temp2=arr[l];
					arr[l]=arr[ir-1];
					arr[ir-1]=temp2;
				}
				if(arr[l-1]>arr[ir-1]){
					temp2=arr[l-1];
					arr[l-1]=arr[ir-1];
					arr[ir-1]=temp2;
				}
				if(arr[l]>arr[l-1]){
					temp2=arr[l];
					arr[l]=arr[l-1];
					arr[l-1]=temp2;
				}
				i=l+1;
				j=ir;
				a=arr[l-1];
				for(;;){
					do
						i++;
					while(arr[i-1]<a);
					do
						j--;
					while(arr[j-1]>a);
					if(j<i)
						break;
					temp2=arr[i-1];
					arr[i-1]=arr[j-1];
					arr[j-1]=temp2;
				}
				arr[l-1]=arr[j-1];
				arr[j-1]=a;
				jstack+=2;
				if(jstack>50)
					System.out.println("NSTACK too small in sort.");
				if(ir-i+1>=j-l){
					istack[jstack-1]=ir;
					istack[jstack-2]=i;
					ir=j-1;
				}else{
					istack[jstack-1]=j-1;
					istack[jstack-2]=l;
					l=i;
				}
			}
		}
	}

	public int[] quicksort_order(float[] arr){
		int i,j,k;
		int ir=arr.length;
		int l=1;
		int jstack=0;
		int[] istack;
		float a,temp2;
		int M=7;
		int b,temp3;

		int[] order=new int[ir];
		for(i=0;i<ir;i++){
			order[i]=i;
		}

		istack=new int[50];
		for(;;){
			if(ir-l<M){
				for(j=l+1;j<=ir;j++){
					a=arr[j-1];
					b=order[j-1];
					for(i=j-1;i>=1;i--){
						if(arr[i-1]<=a)
							break;
						arr[i]=arr[i-1];
						order[i]=order[i-1];
					}
					arr[i]=a;
					order[i]=b;
				}
				if(jstack==0)
					break;
				ir=istack[(jstack--)-1];
				l=istack[(jstack--)-1];
			}else{
				k=(l+ir)>>1;
				temp2=0.0f;
				temp3=0;
				temp2=arr[k-1];
				arr[k-1]=arr[l];
				arr[l]=temp2;
				temp3=order[k-1];
				order[k-1]=order[l];
				order[l]=temp3;
				if(arr[l]>arr[ir-1]){
					temp2=arr[l];
					arr[l]=arr[ir-1];
					arr[ir-1]=temp2;
					temp3=order[l];
					order[l]=order[ir-1];
					order[ir-1]=temp3;
				}
				if(arr[l-1]>arr[ir-1]){
					temp2=arr[l-1];
					arr[l-1]=arr[ir-1];
					arr[ir-1]=temp2;
					temp3=order[l-1];
					order[l-1]=order[ir-1];
					order[ir-1]=temp3;
				}
				if(arr[l]>arr[l-1]){
					temp2=arr[l];
					arr[l]=arr[l-1];
					arr[l-1]=temp2;
					temp3=order[l];
					order[l]=order[l-1];
					order[l-1]=temp3;
				}
				i=l+1;
				j=ir;
				a=arr[l-1];
				b=order[l-1];
				for(;;){
					do
						i++;
					while(arr[i-1]<a);
					do
						j--;
					while(arr[j-1]>a);
					if(j<i)
						break;
					temp2=arr[i-1];
					arr[i-1]=arr[j-1];
					arr[j-1]=temp2;
					temp3=order[i-1];
					order[i-1]=order[j-1];
					order[j-1]=temp3;
				}
				arr[l-1]=arr[j-1];
				order[l-1]=order[j-1];
				arr[j-1]=a;
				order[j-1]=b;
				jstack+=2;
				if(jstack>50)
					System.out.println("NSTACK too small in sort.");
				if(ir-i+1>=j-l){
					istack[jstack-1]=ir;
					istack[jstack-2]=i;
					ir=j-1;
				}else{
					istack[jstack-1]=j-1;
					istack[jstack-2]=l;
					l=i;
				}
			}
		}
		return order;
	}

	public int[] get_order(float[] arr){
		float[] newarr=new float[arr.length];
		System.arraycopy(arr,0,newarr,0,arr.length);
		return quicksort_order(newarr);
	}

	public int[] get_order(int[] arr){
		float[] newarr=new float[arr.length];
		for(int i=0;i<arr.length;i++)
			newarr[i]=arr[i];
		return quicksort_order(newarr);
	}

	public static int[] get_javasort_order(int[] arr){
		List<List<Integer>> list= new ArrayList<List<Integer>>();
		for(int i=0;i<arr.length;i++){
			List<Integer> temp=new ArrayList<Integer>();
			temp.add(Integer.valueOf(arr[i]));
			temp.add(Integer.valueOf(i));
			list.add(temp);
		}
		Collections.sort(list,new Comparator<List<Integer>>(){
			public int compare(List<Integer> o1,List<Integer> o2){
				return o1.get(0).compareTo(o2.get(0));
			}
		});
		int[] order=new int[arr.length];
		for(int i=0;i<arr.length;i++){
			order[i]=list.get(i).get(1).intValue();
		}
		return order;
	}

	public static int[] get_javasort_order(float[] arr){
		List<List<Number>> list=new ArrayList<List<Number>>();
		for(int i=0;i<arr.length;i++){
			List<Number> temp=new ArrayList<Number>();
			temp.add(Float.valueOf(arr[i]));
			temp.add(Integer.valueOf(i));
			list.add(temp);
		}
		Collections.sort(list,new Comparator<List<Number>>(){
			public int compare(List<Number> o1,List<Number> o2){
				Float temp1=(Float)o1.get(0);
				Float temp2=(Float)o2.get(0);
				return temp1.compareTo(temp2);
			}
		});
		int[] order=new int[arr.length];
		for(int i=0;i<arr.length;i++){
			order[i]=list.get(i).get(1).intValue();
		}
		return order;
	}
	
	public static int[] get_javasort_order(String[] arr){
		List<List<Object>> list=new ArrayList<List<Object>>();
		for(int i=0;i<arr.length;i++){
			List<Object> temp=new ArrayList<Object>();
			temp.add(arr[i]);
			temp.add(Integer.valueOf(i));
			list.add(temp);
		}
		Collections.sort(list,new Comparator<List<Object>>(){
			public int compare(List<Object> o1,List<Object> o2){
				String s1=(String)o1.get(0);
				String s2=(String)o2.get(0);
				return s1.compareTo(s2);
			}
		});
		int[] order=new int[arr.length];
		for(int i=0;i<arr.length;i++){
			order[i]=((Integer)list.get(i).get(1)).intValue();
		}
		return order;
	}
	
	public static int[] get_javasort_order(List<String> arr){
		List<List<Object>> list=new ArrayList<List<Object>>();
		for(int i=0;i<arr.size();i++){
			List<Object> temp=new ArrayList<Object>();
			temp.add(arr.get(i));
			temp.add(Integer.valueOf(i));
			list.add(temp);
		}
		Collections.sort(list,new Comparator<List<Object>>(){
			public int compare(List<Object> o1,List<Object> o2){
				String s1=(String)o1.get(0);
				String s2=(String)o2.get(0);
				return s1.compareTo(s2);
			}
		});
		int[] order=new int[arr.size()];
		for(int i=0;i<arr.size();i++){
			order[i]=((Integer)list.get(i).get(1)).intValue();
		}
		return order;
	}

	public static int[] get_javasort_order(long[] arr){
		List<List<Long>> list=new ArrayList<List<Long>>();
		for(int i=0;i<arr.length;i++){
			List<Long> temp=new ArrayList<Long>();
			temp.add(Long.valueOf(arr[i]));
			temp.add(Long.valueOf(i));
			list.add(temp);
		}
		Collections.sort(list,new Comparator<List<Long>>(){
			public int compare(List<Long> o1,List<Long> o2){
				return o1.get(0).compareTo(o2.get(0));
			}
		});
		int[] order=new int[arr.length];
		for(int i=0;i<arr.length;i++){
			order[i]=list.get(i).get(1).intValue();
		}
		return order;
	}

	public static int[] javasort_order(int[] arr){
		float[] newarr=new float[arr.length];
		for(int i=0;i<arr.length;i++)
			newarr[i]=arr[i];
		int[] order=javasort_order(newarr);
		for(int i=0;i<arr.length;i++){
			arr[i]=(int)newarr[i];
		}
		return order;
	}

	public static int[] javasort_order(float[] arr){
		List<List<Number>> list=new ArrayList<List<Number>>();
		for(int i=0;i<arr.length;i++){
			List<Number> temp=new ArrayList<Number>();
			temp.add(Float.valueOf(arr[i]));
			temp.add(Integer.valueOf(i));
			list.add(temp);
		}
		Collections.sort(list,new Comparator<List<Number>>(){
			public int compare(List<Number> o1,List<Number> o2){
				Float temp1=(Float)o1.get(0);
				Float temp2=(Float)o2.get(0);
				return temp1.compareTo(temp2);
			}
		});
		int[] order=new int[arr.length];
		for(int i=0;i<arr.length;i++){
			order[i]=list.get(i).get(1).intValue();
			arr[i]=list.get(i).get(0).floatValue();
		}
		return order;
	}
	
	public static int[] javasort_order(String[] arr){
		List<List<Object>> list=new ArrayList<List<Object>>();
		for(int i=0;i<arr.length;i++){
			List<Object> temp=new ArrayList<Object>();
			temp.add(arr[i]);
			temp.add(Float.valueOf(i));
			list.add(temp);
		}
		Collections.sort(list,new Comparator<List<Object>>(){
			public int compare(List<Object> o1,List<Object> o2){
				String s1=(String)o1.get(0);
				String s2=(String)o2.get(0);
				return s1.compareTo(s2);
			}
		});
		int[] order=new int[arr.length];
		for(int i=0;i<arr.length;i++){
			order[i]=((Float)list.get(i).get(1)).intValue();
			arr[i]=(String)list.get(i).get(0);
		}
		return order;
	}

	public String[] stringsort(String[] list){
		String startval=list[0];
		float[] intlist=new float[list.length];
		for(int i=1;i<list.length;i++){
			intlist[i]=startval.compareTo(list[i]);
		}
		int[] neworder=quicksort_order(intlist);
		String[] temp=new String[list.length];
		for(int i=0;i<list.length;i++){
			temp[i]=list[neworder[i]].substring(0);
		}
		return temp;
	}

	public void sort_string_by(String[] list,float[] sortindex){
		int[] order=get_javasort_order(sortindex);
		String[] temp=new String[list.length];
		for(int i=0;i<list.length;i++){
			temp[i]=list[order[i]].substring(0);
		}
		for(int i=0;i<list.length;i++){
			list[i]=temp[i].substring(0);
		}
	}

	public void sort_string_by(String[] list,long[] sortindex){
		int[] order=get_javasort_order(sortindex);
		String[] temp=new String[list.length];
		for(int i=0;i<list.length;i++){
			temp[i]=list[order[i]].substring(0);
		}
		for(int i=0;i<list.length;i++){
			list[i]=temp[i].substring(0);
		}
	}

	public void sort_object_by(Object[] list,float[] sortindex){
		int[] order=quicksort_order(sortindex);
		Object[] temp=new Object[list.length];
		for(int i=0;i<list.length;i++){
			temp[i]=list[order[i]];
		}
		list=temp;
	}

	public int[] stringsort_order(String[] list,boolean update){
		String startval=list[0];
		float[] intlist=new float[list.length];
		for(int i=1;i<list.length;i++){
			intlist[i]=startval.compareTo(list[i]);
		}
		int[] neworder=quicksort_order(intlist);
		if(update){
			String[] temp=new String[list.length];
			for(int i=0;i<list.length;i++){
				temp[i]=list[neworder[i]].substring(0);
			}
			list=temp;
		}
		return neworder;
	}

}
