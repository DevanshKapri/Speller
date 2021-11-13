// Implements a dictionary's functionality
#include <ctype.h>
#include <stdbool.h>
#include <strings.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "dictionary.h"



// Represents a node in a hash table
typedef struct node
{
    char word[LENGTH + 1];
    struct node *next;
}
node;

// Number of buckets in hash table
const unsigned int N = 10000;
int size_counter = 0;

// Hash table
node *table[N];

// Returns true if word is in dictionary else false
bool check(const char *word)
{
    // Creates copy of word on which hash function can be performed
    int n = strlen(word);
    char word_copy[LENGTH + 1];
    for (int i = 0; i < n; i++)
    {
        word_copy[i] = tolower(word[i]);
    }
    // Adds null terminator to end string
    word_copy[n] = '\0';
    // Initializes index for hashed word
    int h = hash(word_copy);

    node *cursor = table[h];

    while (cursor != NULL)
    {
        if (strcasecmp(cursor->word, word_copy) == 0)
        {
            return true;
        }
        else
        {
            cursor = cursor->next;
        }

    }
    return false;
}

// Hashes word to a number (hash function posted on reddit by delipity)
unsigned int hash(const char *word)
{
    unsigned int hash_output = 0;
    for (int i = 0; word[i] != '\0'; i++)
    {
        hash_output = (hash_output << 2) ^ word[i];
    }
    return hash_output % N;
}

// Loads dictionary into memory, returning true if successful else false
bool load(const char *dictionary)
{
    FILE *file = fopen(dictionary, "r");

    if (file == NULL)
    {
        return false;
    }
    //reading the text word by word and storing the word temporarily before hashing them
    char temp[LENGTH + 1];
    while (fscanf(file, "%s", temp) != EOF)
    {
        //allocated enough memory to store a node
        node *new_node = malloc(sizeof(node));
        //check if malloc has successfully allocated space
        if (new_node == NULL)
        {
            unload();
            fclose(file);
            return false;
        }
        else
        {
            //copies the temp word to the node
            strcpy(new_node->word, temp);
            new_node->next = NULL;

            //hashing the word to give path to the table
            unsigned int hash_index = hash(new_node->word);
            new_node->next = table[hash_index];
            table[hash_index] = new_node;

            size_counter++;
        }
    }
    fclose(file);
    return true;
}

// Returns number of words in dictionary if loaded else 0 if not yet loaded
unsigned int size(void)
{
    return size_counter;
}

// Unloads dictionary from memory, returning true if successful else false
bool unload(void)
{
    for (int i = 0; i < N; i++)
    {
        node *cursor = table[i];
        node *temp = table[i];

        while (cursor != NULL)
        {
            cursor = cursor->next;
            free(temp);
            temp = cursor;
        }
    }
    return true;
}
