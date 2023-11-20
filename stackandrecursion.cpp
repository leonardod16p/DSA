//Solving the Hanoi's Towers Problem recursively with stack

#include <iostream>
#include <vector>
#include <string>

using namespace std;

/*--------------------------------------
 Program Const Areas
 ---------------------------------------*/
const int MAX_SIZE_STACK = 100; // When the stack is defined by an array



/*--------------------------------------
 Data Structure Area
 ---------------------------------------*/

template <class T>
class C_Stack {
    T stack[MAX_SIZE_STACK]; // The Stack is being implemented using a C array 
    int top;
    int id;

public: // INTERFACE

    C_Stack() { top = 0; } // The top is pointing to the first empty position on the Stack
    void push(T data) {stack[top++] = data; }
    T pop(void) { return stack[--top]; }
    bool pilha_esta_vazia(void) { return (top == 0); }
    bool pilha_esta_cheia(void) { return (top == 100); }
    void define_id(int identifier) {id = identifier;}
    int show_id() {return id;}
    int get_top(){
        return top;
    }
    void representation() {
        int i = 0;
        cout << "Tower " << id << ':' << endl;
        cout << '[';
        while (i < top) {
            cout << stack[i];
            i++;
        }
        cout << ']' << endl;
    }
};

//Creating the towers using Stacks

C_Stack<int>* start_scene() {
    C_Stack<int> t1;
    t1.define_id(1);
    C_Stack<int> t2;
    t2.define_id(2);
    C_Stack<int> t3;
    t3.define_id(3);
    C_Stack<int>* towers = new C_Stack<int>[3];
    towers[0] = t1;
    towers[1] = t2;
    towers[2] = t3;
    return towers;
}

C_Stack<int>* towers = start_scene();

//Function resposible for push all the disks on the initial tower 

C_Stack<int> pushing_all_the_disks(C_Stack<int> to_be_pushed, int number_of_disks) {
    int i = number_of_disks;

    while(i > 0){
        to_be_pushed.push(i);
        i = i - 1;
    }
    return to_be_pushed;
}

//Function that take the element on the top of the tower that has start id and puts on the tower with end id

void pop_and_push(int start, int end) {
    int temp = 0;
    cout << "Before: " << endl;
    towers[start].representation();
    towers[end].representation();
    temp = towers[start].pop();
    towers[end].push(temp);
    cout << endl;
    cout << "After: " << endl;
    towers[start].representation();
    towers[end].representation();
    cout << endl;
    cout << endl;

}

//Functions that solves hanoi recursively
void hanoi(int number_of_disks, int start, int end) {
    if (number_of_disks == 1) {
        pop_and_push(start, end);
    }
    else {
        int other_index = 3 - (start + end);
        hanoi(number_of_disks - 1, start, other_index);
        pop_and_push(start, end);
        hanoi(number_of_disks - 1, other_index, end);
    }
}

int main()
{
    int start_at_tower = 0;
    int end_at_tower = 0;
    int number_of_disks = 0;

    cout << "type it the number of disks: ";
    cin >> number_of_disks;
    cout << endl;
    cout << "Type it the identifier of the tower from where you wanna start: ";
    cin >> start_at_tower;
    cout << endl;
    cout << "And the tower you wanna finish: ";
    cin >> end_at_tower;
    cout << endl;
    towers[start_at_tower-1] = pushing_all_the_disks(towers[start_at_tower-1], number_of_disks);


    hanoi(number_of_disks, start_at_tower-1, end_at_tower-1);


    cout << "Final State: " << endl;
    towers[0].representation();
    towers[1].representation();
    towers[2].representation();

    return 0;
}
